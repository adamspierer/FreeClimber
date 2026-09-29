"""
Tests for per-fly tortuosity metrics (scripts/tortuosity.py).

Covers:
  * test_straight_line             -- T==1, S==1, theta==0 for a straight path.
  * test_right_angle_L             -- known T for an L-shape (2 == sqrt(2)/1 ?
                                       no -- L of equal arms: T = 2/sqrt(2) = sqrt(2)).
  * test_closed_loop               -- net_displacement == 0 -> straightness == 0,
                                       tortuosity == NaN (closed bout convention).
  * test_degenerate_short_track    -- <2 points -> all metrics NaN; coincident
                                       points -> path_length == 0, both ratios NaN.
  * test_compute_tortuosity_table  -- end-to-end: synthetic df_filtered with a
                                       'particle' column plus a synthetic df_fng
                                       produce one row per (vial, event, particle).
  * test_cohort_mode_no_particles  -- without a 'particle' column the function
                                       returns an empty, correctly-typed table.

No video decoding required; everything operates on tiny in-memory DataFrames.
"""

import math
import os
import sys

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "scripts"))
import tortuosity as tort  # noqa: E402


# ---------------------------------------------------------------------------
# bout_metrics: pure-function unit tests
# ---------------------------------------------------------------------------
def test_straight_line():
    xs = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    ys = np.zeros_like(xs)
    m = tort.bout_metrics(xs, ys)
    assert m['n_points'] == 5
    assert m['path_length'] == pytest.approx(4.0)
    assert m['net_displacement'] == pytest.approx(4.0)
    assert m['tortuosity'] == pytest.approx(1.0)
    assert m['straightness'] == pytest.approx(1.0)
    assert m['mean_turning_angle_rad'] == pytest.approx(0.0, abs=1e-12)


def test_right_angle_L():
    # Two equal-length arms at 90 degrees: total length = 2, net = sqrt(2).
    xs = np.array([0.0, 1.0, 1.0])
    ys = np.array([0.0, 0.0, 1.0])
    m = tort.bout_metrics(xs, ys)
    assert m['path_length'] == pytest.approx(2.0)
    assert m['net_displacement'] == pytest.approx(math.sqrt(2.0))
    assert m['tortuosity'] == pytest.approx(2.0 / math.sqrt(2.0))
    assert m['straightness'] == pytest.approx(math.sqrt(2.0) / 2.0)
    # Single 90-degree turn -> pi/2.
    assert m['mean_turning_angle_rad'] == pytest.approx(math.pi / 2.0)


def test_closed_loop():
    # Triangle returning to its origin: net displacement == 0.
    xs = np.array([0.0, 1.0, 1.0, 0.0])
    ys = np.array([0.0, 0.0, 1.0, 0.0])
    m = tort.bout_metrics(xs, ys)
    assert m['path_length'] == pytest.approx(2.0 + math.sqrt(2.0))
    assert m['net_displacement'] == pytest.approx(0.0)
    # Convention: tortuosity is NaN for D == 0 (would otherwise be infinite),
    # straightness is exactly 0.
    assert math.isnan(m['tortuosity'])
    assert m['straightness'] == pytest.approx(0.0)
    # Turning angles are still well defined.
    assert not math.isnan(m['mean_turning_angle_rad'])


def test_single_point_track():
    m = tort.bout_metrics([5.0], [5.0])
    assert m['n_points'] == 1
    assert math.isnan(m['path_length'])
    assert math.isnan(m['net_displacement'])
    assert math.isnan(m['tortuosity'])
    assert math.isnan(m['straightness'])
    assert math.isnan(m['mean_turning_angle_rad'])


def test_coincident_points():
    # All points at the same coordinate -> path_length == 0.
    xs = np.array([2.0, 2.0, 2.0, 2.0])
    ys = np.array([3.0, 3.0, 3.0, 3.0])
    m = tort.bout_metrics(xs, ys)
    assert m['n_points'] == 4
    assert m['path_length'] == pytest.approx(0.0)
    assert m['net_displacement'] == pytest.approx(0.0)
    assert math.isnan(m['tortuosity'])
    assert math.isnan(m['straightness'])
    # No non-zero step vectors -> turning angle is undefined.
    assert math.isnan(m['mean_turning_angle_rad'])
    # Vertical efficiency: clamped at 0 when path_length == 0.
    assert m['vertical_efficiency'] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# vertical_efficiency: primary climbing metric (Fix #1)
# ---------------------------------------------------------------------------
def test_vertical_efficiency_perfect_climb():
    """Pure vertical climb -> vertical_efficiency = 1.0."""
    ys = np.linspace(0.0, 10.0, 50)
    xs = np.zeros_like(ys)
    ve = tort.vertical_efficiency(xs, ys)
    assert ve == pytest.approx(1.0, abs=1e-3)


def test_vertical_efficiency_horizontal_walk():
    """Pure horizontal walk -> vertical_efficiency = 0.0."""
    xs = np.linspace(0.0, 10.0, 50)
    ys = np.zeros_like(xs)
    ve = tort.vertical_efficiency(xs, ys)
    assert ve == pytest.approx(0.0, abs=1e-3)


def test_vertical_efficiency_descending():
    """Net descent -> vertical_efficiency clamped at 0.0."""
    xs = np.zeros(50)
    ys = np.linspace(10.0, 0.0, 50)
    ve = tort.vertical_efficiency(xs, ys)
    assert ve == pytest.approx(0.0, abs=1e-3)


def test_vertical_efficiency_diagonal():
    """45-degree diagonal climb -> vertical_efficiency = sin(45 deg) ~ 0.707."""
    xs = np.linspace(0.0, 10.0, 50)
    ys = np.linspace(0.0, 10.0, 50)
    ve = tort.vertical_efficiency(xs, ys)
    assert ve == pytest.approx(math.sin(math.pi / 4.0), abs=1e-3)


# ---------------------------------------------------------------------------
# compute_tortuosity_table: integration unit tests on synthetic DataFrames
# ---------------------------------------------------------------------------
def _make_df_filtered():
    """Synthetic linked detections: two flies in vial 1 over 20 frames.
    Fly 100 climbs straight; fly 101 takes a zig-zag of identical net motion.
    """
    n = 20
    frames = np.arange(n)

    # Fly 100: straight downward y (climbing in image coords).
    fly_a = pd.DataFrame({
        'frame': frames, 't': frames / 25.0, 'vial': 1,
        'x': np.full(n, 100.0),
        'y': np.linspace(0.0, 19.0, n),
        'particle': 100,
    })

    # Fly 101: alternating +/-1 x perturbation, same start/end y.
    x_zig = np.full(n, 110.0) + np.where(frames % 2 == 0, 0.0, 1.0)
    fly_b = pd.DataFrame({
        'frame': frames, 't': frames / 25.0, 'vial': 1,
        'x': x_zig,
        'y': np.linspace(0.0, 19.0, n),
        'particle': 101,
    })

    return pd.concat([fly_a, fly_b], ignore_index=True)


## NOTE: prior FNG-driven bout segmentation tests are intentionally REMOVED
## here per Fix #2 -- bouts are now defined by per-frame vertical velocity,
## independent of FNG events. The replacement tests below exercise that
## algorithm directly.


def test_compute_tortuosity_table_basic_velocity_segmentation():
    """Both flies climb steadily across all 20 frames (~250 mm/s with
    pixel_to_cm=1, frame_rate=25); with velocity_threshold=1 mm/s every
    inter-frame step qualifies -> exactly one bout per particle (2 rows).
    Fly 100's straight path has tortuosity == 1; fly 101's zig-zag > 1."""
    df = _make_df_filtered()
    out = tort.compute_tortuosity_table(
        df, pixel_to_cm=1.0, frame_rate=25.0,
        velocity_threshold=1.0, bout_min_frames=2, bout_min_displacement=0.0,
    )
    assert len(out) == 2
    assert set(out['particle']) == {100, 101}

    straight = out[out['particle'] == 100]
    assert np.allclose(straight['tortuosity'].to_numpy(), 1.0)
    assert np.allclose(straight['straightness'].to_numpy(), 1.0)

    zigzag = out[out['particle'] == 101]
    assert (zigzag['tortuosity'] > 1.0).all()
    assert (zigzag['straightness'] < 1.0).all()


def test_compute_tortuosity_table_no_particle_column():
    """Cohort-mode df_filtered (no 'particle' column) returns an empty table."""
    df = _make_df_filtered().drop(columns=['particle'])
    out = tort.compute_tortuosity_table(df)
    assert out.empty
    assert list(out.columns) == tort.TORTUOSITY_BOUT_COLUMNS


def test_compute_tortuosity_table_missing_required_column_raises():
    df = _make_df_filtered().drop(columns=['x'])
    with pytest.raises(ValueError):
        tort.compute_tortuosity_table(df)


# ---------------------------------------------------------------------------
# Fix #2: velocity-threshold bout segmentation algorithm
# ---------------------------------------------------------------------------
def _piecewise_climb(segments, frame_rate=25.0):
    """Build a synthetic single-particle df_filtered from a list of
    (n_frames, climb_per_frame_px) segments. y starts at 0 and climbs by
    climb_per_frame_px on each frame of the segment; stationary segments use
    climb_per_frame_px == 0.
    """
    rows = []
    f = 0
    y = 0.0
    for n_frames, climb_per_frame in segments:
        for _ in range(n_frames):
            rows.append({'frame': f, 't': f / frame_rate, 'vial': 1,
                         'x': 100.0, 'y': y, 'particle': 1})
            y += climb_per_frame
            f += 1
    return pd.DataFrame(rows)


def test_bout_segmentation_isolates_climbing():
    """Three discrete climbing segments separated by stationary periods
    produce exactly 3 bouts -- segmentation is on velocity, not on FNG."""
    df = _piecewise_climb([
        (15, 1.0),   # climb 1: 15 frames at ~250 mm/s
        (10, 0.0),   # stationary
        (20, 0.8),   # climb 2: 20 frames at ~200 mm/s
        (10, 0.0),   # stationary
        (15, 1.2),   # climb 3: 15 frames at ~300 mm/s
    ])
    out = tort.compute_tortuosity_table(
        df, pixel_to_cm=1.0, frame_rate=25.0,
        velocity_threshold=50.0,   # mm/s: well above 0, below all climbs
        bout_min_frames=2,
        bout_min_displacement=0.0,
    )
    assert len(out) == 3
    assert list(out['bout_idx']) == [1, 2, 3]


def test_bout_segmentation_no_fng_dependency():
    """compute_tortuosity_table operates without any FNG event input. With
    NO df_fng argument it still segments and emits bouts -- confirming
    Fix #2's decoupling from FNG event data."""
    df = _make_df_filtered()
    # Pure positional call: only df_filtered + scalar config. No df_fng.
    out = tort.compute_tortuosity_table(
        df, pixel_to_cm=1.0, frame_rate=25.0,
        velocity_threshold=1.0, bout_min_frames=2, bout_min_displacement=0.0,
    )
    assert not out.empty
    assert len(out) == 2  # one bout per fly
    assert 'bout_idx' in out.columns
    # No FNG-derived 'event_idx' column should appear.
    assert 'event_idx' not in out.columns


def test_short_bouts_filtered():
    """Bouts shorter than bout_min_frames are dropped."""
    df = _piecewise_climb([
        (4, 1.0),    # too short (4 frames)
        (10, 0.0),
        (20, 1.0),   # ok (20 frames)
    ])
    out = tort.compute_tortuosity_table(
        df, pixel_to_cm=1.0, frame_rate=25.0,
        velocity_threshold=10.0,
        bout_min_frames=10,
        bout_min_displacement=0.0,
    )
    assert len(out) == 1
    assert int(out['duration_frames'].iloc[0]) >= 10


def test_smoothing_reduces_noise_inflation(capsys):
    """A perfectly straight vertical climb with Gaussian (sigma=0.5 mm)
    noise added to x and y produces a tortuosity closer to 1.0 with
    Savitzky-Golay smoothing enabled (window=5) than without (window=1)."""
    rng = np.random.default_rng(20260529)
    n = 100
    frames = np.arange(n)
    # Pure climb: y_true = linspace 0..50 px == 0..500 mm (px_to_mm=10 here).
    y_true = np.linspace(0.0, 50.0, n)
    x_true = np.zeros(n)
    # sigma = 0.5 mm in mm-space == 0.05 px (with pixel_to_cm=1.0 -> 1 px == 10 mm)
    sigma_px = 0.05
    xs = x_true + rng.normal(0.0, sigma_px, n)
    ys = y_true + rng.normal(0.0, sigma_px, n)

    df = pd.DataFrame({'frame': frames, 'vial': 1,
                       'x': xs, 'y': ys, 'particle': 1})

    # No smoothing: window = 1 disables Sav-Gol per _normalize_savgol_window.
    out_raw = tort.compute_tortuosity_table(
        df, pixel_to_cm=1.0, frame_rate=25.0,
        velocity_threshold=0.0,
        bout_min_frames=2, bout_min_displacement=0.0,
        smoothing_window=1,
    )
    out_smoothed = tort.compute_tortuosity_table(
        df, pixel_to_cm=1.0, frame_rate=25.0,
        velocity_threshold=0.0,
        bout_min_frames=2, bout_min_displacement=0.0,
        smoothing_window=5,
    )

    assert not out_raw.empty and not out_smoothed.empty
    t_raw = float(out_raw['tortuosity'].iloc[0])
    t_smoothed = float(out_smoothed['tortuosity'].iloc[0])

    # Smoothing should pull T strictly closer to the ideal value of 1.0.
    assert abs(t_smoothed - 1.0) < abs(t_raw - 1.0), (
        'smoothing did not reduce noise-driven tortuosity inflation: '
        'raw T=%.6f smoothed T=%.6f' % (t_raw, t_smoothed)
    )


def test_smoothing_skipped_for_short_trajectories(capsys):
    """A particle with 3 frames and smoothing_window=5 produces a clear
    stdout warning and computation proceeds using the raw coordinates."""
    df = pd.DataFrame({
        'frame': [0, 1, 2], 'vial': 1,
        'x': [0.0, 0.0, 0.0], 'y': [0.0, 1.0, 2.0],
        'particle': 42,
    })
    out = tort.compute_tortuosity_table(
        df, pixel_to_cm=1.0, frame_rate=25.0,
        velocity_threshold=0.0,
        bout_min_frames=2, bout_min_displacement=0.0,
        smoothing_window=5,
    )
    captured = capsys.readouterr().out
    assert 'skipping Savitzky-Golay smoothing' in captured, (
        'expected a warning about skipped smoothing; got:\n%s' % captured)
    assert 'particle 42' in captured
    # Computation still produced a bout (raw coords are a straight climb).
    assert len(out) == 1
    assert int(out['particle'].iloc[0]) == 42


def test_low_displacement_bouts_filtered():
    """Bouts whose net vertical displacement is below the configured
    threshold (mm) are dropped."""
    # Two climbs: one rises ~1 mm total, one rises ~10 mm total.
    df = _piecewise_climb([
        (15, 0.005),   # ~15 * 0.005 px = 0.075 px; in mm with pixel_to_cm=1
                       # that is 0.75 mm -- below the 5 mm threshold.
        (10, 0.0),
        (15, 0.10),    # ~15 * 0.10 px = 1.5 px = 15.0 mm -- above 5 mm.
    ])
    out = tort.compute_tortuosity_table(
        df, pixel_to_cm=1.0, frame_rate=25.0,
        # Threshold low enough that the slow climb is also "climbing"; the
        # FILTER then drops the low-displacement bout.
        velocity_threshold=0.5,
        bout_min_frames=2,
        bout_min_displacement=5.0,   # mm
    )
    assert len(out) == 1
    assert float(out['vertical_displacement_mm'].iloc[0]) >= 5.0


# ---------------------------------------------------------------------------
# Integration with the detector hook -- writes <video>.tortuosity.csv
# ---------------------------------------------------------------------------
import types  # noqa: E402

sys.path.insert(0, os.path.join(PROJECT_ROOT, "scripts"))
import detector_fng as dfng  # noqa: E402


def _bind(det, *names):
    # Always bind the output relabeling helpers: every CSV-writing method now
    # passes its frame through _relabel_vial_col (which calls _vial_label).
    for name in tuple(names) + ("_vial_label", "_relabel_vial_col"):
        setattr(det, name, types.MethodType(getattr(dfng.detector, name), det))
    return det


def test_detector_compute_tortuosity_writes_csv(tmp_path):
    """detector.compute_tortuosity() writes BOTH *.tortuosity_bouts.csv and
    *.tortuosity_particle.csv next to the other outputs. The detector is
    invoked WITHOUT df_fng to exercise Fix #2's FNG-independence."""
    df = _make_df_filtered()

    det = types.SimpleNamespace(
        debug=False,
        df_filtered=df,
        # df_fng INTENTIONALLY ABSENT -- compute_tortuosity must not read it.
        name_nosuffix=os.path.join(str(tmp_path), "synthetic"),
        pixel_to_cm=1.0,
        frame_rate=25.0,
        tortuosity_velocity_threshold=1.0,
        tortuosity_bout_min_frames=2,
        tortuosity_bout_min_displacement=0.0,
    )
    _bind(det, "compute_tortuosity")
    det.compute_tortuosity()

    bouts_path = os.path.join(str(tmp_path), "synthetic.tortuosity_bouts.csv")
    particle_path = os.path.join(str(tmp_path), "synthetic.tortuosity_particle.csv")
    assert os.path.exists(bouts_path)
    assert os.path.exists(particle_path)

    bouts = pd.read_csv(bouts_path)
    assert list(bouts.columns) == tort.TORTUOSITY_BOUT_COLUMNS
    assert len(bouts) == 2   # one bout per fly under uniform climbing
    assert 'vertical_efficiency' in bouts.columns

    particle = pd.read_csv(particle_path)
    assert list(particle.columns) == tort.TORTUOSITY_PARTICLE_COLUMNS
    assert 'median_vertical_efficiency' in particle.columns
    assert len(particle) == 2  # two particles -> two rows


def test_detector_compute_tortuosity_cohort_mode_writes_empty_csv(tmp_path):
    """Cohort mode (no 'particle' column) still emits header-only csvs so
    downstream tooling doesn't have to special-case missing files."""
    df = _make_df_filtered().drop(columns=['particle'])
    det = types.SimpleNamespace(
        debug=False,
        df_filtered=df,
        name_nosuffix=os.path.join(str(tmp_path), "cohort"),
        pixel_to_cm=1.0, frame_rate=25.0,
    )
    _bind(det, "compute_tortuosity")
    det.compute_tortuosity()

    bouts_path = os.path.join(str(tmp_path), "cohort.tortuosity_bouts.csv")
    particle_path = os.path.join(str(tmp_path), "cohort.tortuosity_particle.csv")
    assert os.path.exists(bouts_path)
    assert os.path.exists(particle_path)

    bouts = pd.read_csv(bouts_path)
    assert bouts.empty
    assert list(bouts.columns) == tort.TORTUOSITY_BOUT_COLUMNS
    particle = pd.read_csv(particle_path)
    assert particle.empty
    assert list(particle.columns) == tort.TORTUOSITY_PARTICLE_COLUMNS
