#!/usr/bin/env python
# -*- coding: utf-8 -*-

## File name : tortuosity.py
## Purpose   : Per-fly trajectory tortuosity metrics for FreeClimber-FNG.
##             Computed from linked tracks (analysis_mode == 'individual').

"""
Per-fly tortuosity metrics for FreeClimber-FNG.

Four metrics are computed per climbing bout, per fly:

  1. tortuosity (T)            classical tortuosity = path_length / net_displacement
                               T >= 1 for any non-zero net displacement; larger
                               values mean a more meandering path.
  2. straightness (S)          straightness index   = net_displacement / path_length
                               S in [0, 1]; S = 1 is a perfectly straight path.
                               (S == 1 / T; we record both because each is
                                conventional in different sub-fields and the
                                inverse is undefined when path_length == 0.)
  3. vertical_efficiency       max(0, y_end - y_start) / path_length, in [0, 1].
                               1.0 == every unit of path translates into one unit
                               of upward displacement; 0.0 == no net climb (or a
                               net descent). The PRIMARY climbing metric.
  4. mean_turning_angle (rad)  AUXILIARY path-shape metric: mean absolute
                               turning angle between consecutive step vectors,
                               in radians. Kept for closed-loop bouts and as a
                               directional-persistence indicator.

BOUT SEGMENTATION  (Fix #2 -- INDEPENDENT of FNG events)
  For each (vial, particle) trajectory, bouts are detected by per-frame
  vertical climbing velocity. A bout is a maximal run of consecutive frames
  whose step-wise vertical velocity exceeds velocity_threshold (mm/s).
  Bouts are terminated when the velocity falls to or below the threshold for
  any single inter-frame step. Bouts shorter than bout_min_frames frames OR
  with a net vertical displacement below bout_min_displacement mm are then
  dropped.

  compute_tortuosity_table does NOT consult FNG event data; it has no
  dependency on _detect_fng_series, self.fng_events, or any FNG output.

Y-CONVENTION
  detector.step_5 inverts y via detector.invert_y (detector_fng.py:1492), so
  in self.df_filtered "climbing" means INCREASING y. Confirmed empirically
  against the cohort climbing-slope sign reported by
  detector.local_linear_regression on the FNG synthetic validation fixture
  (slope is positive, +4.7428 px/frame, during the climb; mean y rises from
  ~7.8 in early frames to ~155 in late frames).

UNITS
  vertical_efficiency, tortuosity and straightness are unit-free, but are
  defined for path_length and (y_end - y_start) measured in the SAME linear
  units. compute_tortuosity in detector_fng.py converts pixel coordinates to
  MILLIMETRES (using self.pixel_to_cm, which is pixels-per-cm) before
  invoking these helpers, so reported path lengths and net displacements are
  in mm and velocity thresholds are in mm/s.
"""

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

TORTUOSITY_BOUT_COLUMNS = [
    'vial', 'particle', 'bout_idx',
    'frame_start', 'frame_end', 'n_points',
    'duration_frames', 'duration_sec',
    'path_length_mm', 'net_displacement_mm',
    'vertical_displacement_mm',
    'tortuosity', 'straightness', 'vertical_efficiency',
    'mean_turning_angle_rad',
]

# Backward-compat alias for any external callers / tests that imported
# TORTUOSITY_COLUMNS from earlier versions of this module.
TORTUOSITY_COLUMNS = TORTUOSITY_BOUT_COLUMNS

TORTUOSITY_PARTICLE_COLUMNS = [
    'vial', 'particle', 'n_bouts', 'median_vertical_efficiency',
]


# ---------------------------------------------------------------------------
# Pure metric helpers (per-bout)
# ---------------------------------------------------------------------------
def path_length(xs, ys):
    """Cumulative Euclidean distance along the (xs, ys) polyline.

    Returns 0.0 for sequences with fewer than 2 points.
    """
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    if xs.size < 2:
        return 0.0
    dx = np.diff(xs)
    dy = np.diff(ys)
    return float(np.sum(np.sqrt(dx * dx + dy * dy)))


def net_displacement(xs, ys):
    """Straight-line distance from first to last point of the (xs, ys) polyline.

    Returns 0.0 for sequences with fewer than 2 points.
    """
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    if xs.size < 2:
        return 0.0
    return float(np.hypot(xs[-1] - xs[0], ys[-1] - ys[0]))


def vertical_efficiency(xs, ys):
    """Fraction of (smoothed) path length that maps to net upward climb.

    Definition (per spec): vertical_efficiency = max(0, y_end - y_start) / L
    where L is the smoothed path length over (xs, ys) and (xs, ys) are in
    millimetres. Pure function: callers in compute_tortuosity smooth the
    coordinates (Savitzky-Golay, polyorder=2) and convert to mm before
    invoking this helper.

    Y-convention: in self.df_filtered, "climbing" means INCREASING y -- see
    detector.invert_y at detector_fng.py:745-757 and the call site at
    detector_fng.py:1492. Empirically confirmed against the sign of the
    cohort climbing slope from detector.local_linear_regression on the
    synthetic FNG validation fixture (positive slope during climb).

    Range: [0, 1]. Returns 0.0 when path_length == 0 or when y_end <= y_start.
    """
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    if xs.size < 2:
        return 0.0
    L = path_length(xs, ys)
    if L == 0.0:
        return 0.0
    rise = float(ys[-1] - ys[0])
    if rise <= 0.0:
        return 0.0
    eff = rise / L
    # Numerical guard: rise can never exceed total path length, but clip
    # to keep the contract [0, 1] watertight in the face of FP noise.
    if eff > 1.0:
        return 1.0
    return eff


def mean_absolute_turning_angle(xs, ys):
    """Mean absolute turning angle (radians) between consecutive step vectors.

    For a polyline with n points there are n-1 steps and n-2 turning angles.
    Steps with zero length are skipped (they have no defined direction).
    Returns NaN if fewer than 2 well-defined turning angles can be computed.

    Auxiliary path-shape metric: not the primary climbing metric (see
    vertical_efficiency) but kept for closed-loop bouts where T and S
    degenerate and as a directional-persistence indicator.
    """
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    if xs.size < 3:
        return float('nan')

    vx = np.diff(xs)
    vy = np.diff(ys)
    norms = np.sqrt(vx * vx + vy * vy)

    # Drop zero-length steps; a zero-length step has no defined direction.
    valid = norms > 0
    if valid.sum() < 2:
        return float('nan')
    vx = vx[valid]
    vy = vy[valid]
    norms = norms[valid]

    # Cosine of turning angle between consecutive step vectors.
    cos_theta = (vx[:-1] * vx[1:] + vy[:-1] * vy[1:]) / (norms[:-1] * norms[1:])
    # Floating-point can push |cos| just past 1, which would NaN arccos.
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    angles = np.arccos(cos_theta)
    return float(np.mean(np.abs(angles)))


def bout_metrics(xs, ys):
    """Compute the per-bout tortuosity metric bundle for a single polyline.

    Inputs:
      xs, ys : 1-D arrays of equal length, ordered by frame (ascending). For
               spec conformance compute_tortuosity passes smoothed (mm)
               coordinates; this pure function does not require it.

    Returns:
      dict with keys: n_points, path_length, net_displacement, tortuosity,
                     straightness, vertical_efficiency, mean_turning_angle_rad.

    Degenerate handling:
      * fewer than 2 points -> all metrics NaN (vertical_efficiency = 0).
      * net_displacement == 0 (closed loop) -> tortuosity is NaN
        (path_length / 0 has no useful value); straightness is 0.
      * path_length == 0 (all points coincide) -> straightness is NaN
        and tortuosity is NaN; vertical_efficiency is 0.
    """
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    n = int(xs.size)

    if n < 2:
        return {
            'n_points': n,
            'path_length': float('nan'),
            'net_displacement': float('nan'),
            'tortuosity': float('nan'),
            'straightness': float('nan'),
            'vertical_efficiency': 0.0,
            'mean_turning_angle_rad': float('nan'),
        }

    L = path_length(xs, ys)
    D = net_displacement(xs, ys)

    if L == 0.0:
        T = float('nan')
        S = float('nan')
    elif D == 0.0:
        # Closed loop: classical tortuosity diverges; straightness is exactly 0.
        T = float('nan')
        S = 0.0
    else:
        T = L / D
        S = D / L

    theta = mean_absolute_turning_angle(xs, ys)
    ve = vertical_efficiency(xs, ys)

    return {
        'n_points': n,
        'path_length': L,
        'net_displacement': D,
        'tortuosity': T,
        'straightness': S,
        'vertical_efficiency': ve,
        'mean_turning_angle_rad': theta,
    }


# ---------------------------------------------------------------------------
# Trajectory smoothing (Fix #4)
# ---------------------------------------------------------------------------
def _normalize_savgol_window(window):
    """Coerce a user-supplied window length to a positive odd integer.

    scipy.signal.savgol_filter requires window_length to be odd and at least
    polyorder + 1. Returns the normalized window length; returns 0 to signal
    "no smoothing" when window <= 1.
    """
    try:
        w = int(window)
    except (TypeError, ValueError):
        return 0
    if w <= 1:
        return 0
    if w % 2 == 0:
        w += 1
    return w


def smooth_xy(xs, ys, window, particle_label=None):
    """Apply Savitzky-Golay smoothing (polyorder=2, mode='nearest') to x and y.

    If the input has fewer frames than the (normalized) window length the
    smoothing is SKIPPED for this particle, a warning line is printed to
    stdout, and the raw arrays are returned unchanged (per spec).

    Returns (xs_smoothed, ys_smoothed, smoothing_applied) where the third
    element is a bool indicating whether smoothing actually ran.
    """
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    w = _normalize_savgol_window(window)

    if w == 0:
        return xs, ys, False

    n = int(xs.size)
    if n < w:
        label = '' if particle_label is None else ' for particle %s' % particle_label
        print('   !! tortuosity: skipping Savitzky-Golay smoothing%s '
              '(have %d frames, need %d); using raw coordinates'
              % (label, n, w))
        return xs, ys, False

    polyorder = min(2, w - 1)
    xs_s = savgol_filter(xs, window_length=w, polyorder=polyorder, mode='nearest')
    ys_s = savgol_filter(ys, window_length=w, polyorder=polyorder, mode='nearest')
    return xs_s, ys_s, True


# ---------------------------------------------------------------------------
# Bout segmentation (Fix #2: velocity-threshold, NOT FNG-driven)
# ---------------------------------------------------------------------------
def _segment_climbing_bouts(frames, ys_mm, frame_rate, velocity_threshold):
    """Identify climbing bouts in a single trajectory by vertical velocity.

    A bout is a maximal contiguous run of consecutive (frame_k, frame_{k+1})
    steps whose step-wise vertical velocity exceeds velocity_threshold mm/s.
    Note: contiguous in *position within the per-particle frame sequence*,
    NOT in absolute frame index -- so a single skipped detection terminates
    a bout, which is the conservative behavior for noisy/intermittent
    detections.

    Inputs:
      frames             : 1-D int array, sorted ascending, the per-frame
                           indices for one particle.
      ys_mm              : 1-D float array, same length, vertical positions
                           in millimetres (climbing = increasing).
      frame_rate         : frames per second.
      velocity_threshold : minimum vertical step velocity in mm/s.

    Returns:
      list of (start_pos, end_pos) inclusive tuples giving the index range
      in frames/ys_mm covered by each bout. A bout spanning steps k..m
      covers positions k..m+1 (m-k+2 frames).
    """
    n = int(len(frames))
    if n < 2:
        return []

    fr = float(frame_rate) if frame_rate else 1.0
    df_frames = np.diff(frames).astype(float)
    dy = np.diff(ys_mm).astype(float)

    # Step duration in seconds. Guard against zero-duration steps (duplicate
    # frame numbers should never occur, but if they do treat the step as
    # below threshold).
    dt_sec = np.where(df_frames > 0, df_frames / fr, np.inf)
    vy = np.where(np.isfinite(dt_sec), dy / dt_sec, 0.0)
    climbing_step = vy > float(velocity_threshold)

    bouts = []
    in_bout = False
    bout_start = 0
    for k, climbing in enumerate(climbing_step):
        if climbing:
            if not in_bout:
                in_bout = True
                bout_start = k
        else:
            if in_bout:
                bouts.append((bout_start, k))  # last climbing step was k-1; covers frames bout_start..k
                in_bout = False
    if in_bout:
        bouts.append((bout_start, len(climbing_step)))
    return bouts


def compute_tortuosity_table(df_filtered,
                             pixel_to_cm=1.0,
                             frame_rate=1.0,
                             velocity_threshold=0.0,
                             bout_min_frames=1,
                             bout_min_displacement=0.0,
                             smoothing_window=1):
    """Build the per-(vial, particle, bout) tortuosity-bout table.

    Pipeline per particle (Fix #4):
      1. Convert (x, y) from pixels to millimetres using pixel_to_cm.
      2. Apply Savitzky-Golay smoothing (polyorder=2, mode='nearest') with
         the configured smoothing_window. Particles with fewer frames than
         the window are skipped (warning logged) and use raw coordinates.
      3. Detect climbing bouts from the SMOOTHED y-velocity.
      4. Drop bouts below bout_min_frames or bout_min_displacement.
      5. Compute per-bout path_length, net_displacement,
         vertical_displacement, vertical_efficiency, tortuosity,
         straightness, mean_turning_angle from the SMOOTHED coordinates.

    Bout segmentation is INDEPENDENT of FNG events (Fix #2).

    Inputs:
      df_filtered           : DataFrame with columns frame, vial, x, y,
                              particle. (x, y) are in pixels with the
                              df_filtered y-inversion convention (climbing =
                              increasing y). Cohort-mode input (no 'particle'
                              column) returns an empty table.
      pixel_to_cm           : conversion factor (pixels per cm). Converts
                              x, y from pixels to millimetres internally.
      frame_rate            : frames per second.
      velocity_threshold    : minimum vertical step velocity, in mm/s, for
                              a step to count as a climbing step.
      bout_min_frames       : drop bouts shorter than this many frames.
      bout_min_displacement : drop bouts whose net vertical displacement is
                              less than this many millimetres.
      smoothing_window      : Savitzky-Golay window length in frames
                              (polyorder=2). <=1 disables smoothing; even
                              values are rounded up to the next odd.

    Returns:
      DataFrame with columns TORTUOSITY_BOUT_COLUMNS. Empty (correctly-typed
      headers only) when there is nothing to compute.
    """
    if df_filtered is None or df_filtered.empty or 'particle' not in df_filtered.columns:
        return pd.DataFrame(columns=TORTUOSITY_BOUT_COLUMNS)

    needed = {'frame', 'vial', 'x', 'y', 'particle'}
    missing = needed - set(df_filtered.columns)
    if missing:
        raise ValueError('df_filtered missing required column(s): %s' % sorted(missing))

    p2cm = float(pixel_to_cm) if pixel_to_cm else 1.0
    px_to_mm = 10.0 / p2cm  # cm * 10 = mm, and pixels / (px per cm) = cm

    records = []

    for (vial, particle), dfp in df_filtered.groupby(['vial', 'particle']):
        dfp = dfp.sort_values('frame').reset_index(drop=True)
        frames = dfp['frame'].to_numpy()
        xs_mm = dfp['x'].to_numpy(dtype=float) * px_to_mm
        ys_mm = dfp['y'].to_numpy(dtype=float) * px_to_mm

        # Savitzky-Golay smoothing of (x, y) per particle, in mm. ALL
        # downstream computations -- velocity, bout segmentation, path
        # length, net displacement, vertical_efficiency, tortuosity,
        # straightness, turning angle -- use the smoothed arrays.
        xs_mm, ys_mm, _ = smooth_xy(xs_mm, ys_mm, smoothing_window,
                                    particle_label=int(particle))

        bouts = _segment_climbing_bouts(frames, ys_mm,
                                        frame_rate=frame_rate,
                                        velocity_threshold=velocity_threshold)
        bout_ord = 0
        for start_pos, end_pos in bouts:
            # end_pos is the exclusive position of the last *climbing step*;
            # the bout covers frame positions start_pos .. end_pos (inclusive),
            # i.e. end_pos - start_pos + 1 points.
            stop_pos = end_pos
            xb = xs_mm[start_pos:stop_pos + 1]
            yb = ys_mm[start_pos:stop_pos + 1]
            fb = frames[start_pos:stop_pos + 1]

            if len(fb) < 2:
                continue

            duration_frames = int(fb[-1] - fb[0] + 1)
            vertical_displacement_mm = float(yb[-1] - yb[0])

            if duration_frames < int(bout_min_frames):
                continue
            if vertical_displacement_mm < float(bout_min_displacement):
                continue

            m = bout_metrics(xb, yb)
            bout_ord += 1
            records.append({
                'vial': int(vial),
                'particle': int(particle),
                'bout_idx': int(bout_ord),
                'frame_start': int(fb[0]),
                'frame_end': int(fb[-1]),
                'n_points': int(len(fb)),
                'duration_frames': duration_frames,
                'duration_sec': round(duration_frames / float(frame_rate), 6) if frame_rate else float('nan'),
                'path_length_mm': m['path_length'],
                'net_displacement_mm': m['net_displacement'],
                'vertical_displacement_mm': vertical_displacement_mm,
                'tortuosity': m['tortuosity'],
                'straightness': m['straightness'],
                'vertical_efficiency': m['vertical_efficiency'],
                'mean_turning_angle_rad': m['mean_turning_angle_rad'],
            })

    if not records:
        return pd.DataFrame(columns=TORTUOSITY_BOUT_COLUMNS)
    return pd.DataFrame.from_records(records, columns=TORTUOSITY_BOUT_COLUMNS)


def compute_particle_table(bouts):
    """Build the per-particle aggregate table from a per-bout table.

    Inputs:
      bouts : DataFrame with the TORTUOSITY_BOUT_COLUMNS schema.

    Returns:
      DataFrame with columns TORTUOSITY_PARTICLE_COLUMNS:
        vial, particle, n_bouts, median_vertical_efficiency.

      The spec for the particle CSV currently names median_vertical_efficiency
      explicitly; the file's primary purpose is per-fly aggregation of the
      primary climbing metric.
    """
    if bouts is None or bouts.empty:
        return pd.DataFrame(columns=TORTUOSITY_PARTICLE_COLUMNS)

    grouped = (bouts
               .groupby(['vial', 'particle'], as_index=False)
               .agg(n_bouts=('vial', 'size'),
                    median_vertical_efficiency=('vertical_efficiency', 'median')))
    return grouped[TORTUOSITY_PARTICLE_COLUMNS]
