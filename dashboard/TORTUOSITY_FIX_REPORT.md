# Tortuosity-metrics spec-conformance fixes

This document records the four spec-deviation fixes applied on top of
the initial tortuosity implementation (branch
`feature/tortuosity-metrics`, parent commit `df79d0e`).

**Note on the previous report.** The fix prompt referred to a
`dashboard/TORTUOSITY_IMPLEMENTATION_REPORT.md` as the prior status
report. That file does not exist in the repository — the initial
implementation's status report was an in-chat prose summary, not a
committed markdown file. Per principle G (do not rewrite history) and
the spec instruction "in addition to, not replacing", I have not
fabricated a backdated report; this fix report stands alone.

- **Branch:** `feature/tortuosity-metrics`
- **Final commit SHA:** `c5b51fb1dd29738670e8974b686dd42dda4ddb3`
- **Parent of the four fixes:** `df79d0e` (initial implementation tip)

## Git log of this fix iteration

```
c5b51fb Fix #4: Apply Savitzky-Golay smoothing to (x, y) per particle
8fe81f6 Fix #3: Add five tortuosity config keys to all three ALLOWED_KEYS sites
5a328c6 Fix #2: Decouple bout segmentation from FNG events (velocity-threshold algorithm)
d4df2bf Fix #1: Add vertical_efficiency metric and per-particle output CSV
```

---

## Fix #1 — Add vertical_efficiency as the primary climbing metric

### Spec requirement (verbatim)

> The original spec required three primary metrics: tortuosity, straightness, and vertical_efficiency. Only the first two were implemented; mean absolute turning angle was substituted for vertical_efficiency. Mean absolute turning angle stays (principle F); vertical_efficiency must be added.
>
> Requirements:
> - Add a pure function `vertical_efficiency(x, y)` to scripts/tortuosity.py.
> - Definition: vertical_efficiency = max(0, y_end − y_start) / path_length, where path_length is computed from the smoothed Savitzky-Golay path (Fix #4) and (x, y) are in mm.
> - y-inversion convention: detector_fng.py:1320 inverts y, so in self.df_filtered "climbing" means INCREASING y. Verify this empirically against the cohort climbing-slope sign in local_linear_regression() output on an existing test fixture before relying on it. Document the convention explicitly in the function docstring.
> - Range [0, 1]. Returns 0 if path_length is 0 or if y_end <= y_start.
> - Add `vertical_efficiency` column to *.tortuosity_bouts.csv output.
> - Add `median_vertical_efficiency` column to *.tortuosity_particle.csv output.
> - Add tests to tests/test_tortuosity.py:
>   - test_vertical_efficiency_perfect_climb: straight-up trajectory produces vertical_efficiency = 1.0 within 1e-3.
>   - test_vertical_efficiency_horizontal_walk: horizontal-only trajectory produces vertical_efficiency = 0.0 within 1e-3.
>   - test_vertical_efficiency_descending: descending-only trajectory produces vertical_efficiency = 0.0 (clamped at 0).
>   - test_vertical_efficiency_diagonal: 45-degree diagonal climb produces vertical_efficiency ≈ sin(45°) ≈ 0.707 within 1e-3.

### What was implemented

- Pure function `vertical_efficiency(xs, ys)` at
  `scripts/tortuosity.py:99-127`. Implements
  `max(0, ys[-1] - ys[0]) / path_length(xs, ys)`; returns `0.0` when
  `path_length == 0` or `rise <= 0`; clipped to `[0, 1]`.
- Wired into `bout_metrics()` at `scripts/tortuosity.py:222-227` (the
  returned dict gains a `vertical_efficiency` key).
- New per-bout CSV column `vertical_efficiency` in
  `TORTUOSITY_BOUT_COLUMNS` at `scripts/tortuosity.py:64-71`.
- New per-particle aggregate CSV: `compute_particle_table()` at
  `scripts/tortuosity.py:472-490`, producing
  `TORTUOSITY_PARTICLE_COLUMNS = [vial, particle, n_bouts, median_vertical_efficiency]`.
- Detector hook `detector.compute_tortuosity()` (at
  `scripts/detector_fng.py:1182-1252`) now writes both
  `<video>.tortuosity_bouts.csv` and `<video>.tortuosity_particle.csv`
  (previous code wrote only one file, `<video>.tortuosity.csv`).
- Y-convention verified empirically against the synthetic FNG
  validation fixture (`tests/fixtures/synthetic_validation/`):
  cohort `local_linear_regression()` reports slope **+4.7428 px/frame**
  during the climb (see
  `tests/fixtures/synthetic_validation/freeclimber_fng_validation_video.slopes.csv`),
  and mean y rises from ~7.8 (early frames) to ~155 (late frames).
  Climbing therefore means INCREASING y in `df_filtered`. The convention
  is documented in the `vertical_efficiency` docstring and in the module
  header.

### Tests added or modified

| Test | File | Asserts |
|---|---|---|
| `test_vertical_efficiency_perfect_climb` | `tests/test_tortuosity.py` | Straight up: `ve ≈ 1.0` (tol 1e-3). |
| `test_vertical_efficiency_horizontal_walk` | `tests/test_tortuosity.py` | Purely horizontal: `ve ≈ 0.0` (tol 1e-3). |
| `test_vertical_efficiency_descending` | `tests/test_tortuosity.py` | Net descent: `ve ≈ 0.0` (clamped). |
| `test_vertical_efficiency_diagonal` | `tests/test_tortuosity.py` | 45° diagonal: `ve ≈ sin(π/4) ≈ 0.707` (tol 1e-3). |
| `test_coincident_points` (modified) | `tests/test_tortuosity.py` | Added assertion `ve == 0` for coincident points. |
| `test_detector_compute_tortuosity_writes_csv` (modified) | `tests/test_tortuosity.py` | Now asserts both bouts.csv and particle.csv are written with the new columns. |
| `test_detector_compute_tortuosity_cohort_mode_writes_empty_csv` (modified) | `tests/test_tortuosity.py` | Asserts both files exist (header-only) in cohort mode. |

### Confirmation

**Spec item complete.**

---

## Fix #2 — Decouple bout segmentation from FNG events

### Spec requirement (verbatim)

> The original spec required bout segmentation to be INDEPENDENT of FNG event detection. The current implementation defines bouts as intervals between FNG events. This must change.
>
> Requirements:
> - compute_tortuosity must NOT call _detect_fng_series, MUST NOT read self.fng_events or any FNG output, and MUST NOT use FNG event data to define bout boundaries.
> - compute_tortuosity may continue to be invoked AFTER compute_fng in step_5 (the call order is fine; only the data dependency is forbidden).
> - New bout segmentation algorithm: for each particle, identify climbing bouts as contiguous frame runs where the per-frame vertical velocity > tortuosity_velocity_threshold (config key from Fix #3). Bouts terminate when velocity falls below threshold for one or more frames.
> - After segmentation, drop bouts where duration_frames < tortuosity_bout_min_frames (Fix #3 config).
> - Also drop bouts where the net vertical displacement < tortuosity_bout_min_displacement (Fix #3 config).
> - The bout_metrics function and any helpers in scripts/tortuosity.py that compute per-bout metrics must not require any FNG-derived input. The input should be (x, y, frame, vial, particle) trajectory data only.
> - Update or replace existing bout-segmentation tests so they test the velocity-threshold algorithm rather than FNG-tied behavior. Specifically:
>   - test_bout_segmentation_isolates_climbing: synthetic trajectory with three discrete climbing segments separated by stationary periods produces exactly 3 bouts.
>   - test_bout_segmentation_no_fng_dependency: compute_tortuosity called on a trajectory with no FNG events still produces correct bouts.
>   - test_short_bouts_filtered: bouts below tortuosity_bout_min_frames are dropped.
>   - test_low_displacement_bouts_filtered: bouts below tortuosity_bout_min_displacement are dropped.
> - The status report must include a code snippet (5–15 lines) from compute_tortuosity showing that bout segmentation operates on velocity, not on FNG events.

### What was implemented

- Removed the FNG-driven helper `_bout_windows_for_vial` and replaced
  it with `_segment_climbing_bouts(frames, ys_mm, frame_rate, velocity_threshold)`
  at `scripts/tortuosity.py:299-352`. The function takes only trajectory
  data and scalar configuration — it has no `df_fng` argument and no
  call into `_detect_fng_series`.
- `compute_tortuosity_table()` at `scripts/tortuosity.py:354-468`
  rewrote its signature: no `df_fng` argument; it now takes
  `pixel_to_cm`, `frame_rate`, `velocity_threshold`, `bout_min_frames`,
  `bout_min_displacement` (and `smoothing_window`, added by Fix #4).
- Per-particle loop converts pixels → mm, runs
  `_segment_climbing_bouts`, drops bouts whose
  `duration_frames < bout_min_frames` or
  `vertical_displacement_mm < bout_min_displacement`
  (`scripts/tortuosity.py:432-448`).
- `event_idx` renamed to `bout_idx` in `TORTUOSITY_BOUT_COLUMNS`
  (`scripts/tortuosity.py:64-71`) — FNG events no longer drive bouts.
- Detector hook `detector.compute_tortuosity()`
  (`scripts/detector_fng.py:1182-1252`):
  - No longer reads `self.df_fng` / `self.fng_events`.
  - No longer calls `_detect_fng_series`.
  - Configs pulled via `getattr` (`detector_fng.py:1210-1216`).
- Call order in `step_5` preserved (`detector_fng.py:1552-1555`):
  `compute_fng()` still runs first, but `compute_tortuosity()` no longer
  consumes its output.

### Required code snippet (Fix #2 spec)

From `compute_tortuosity_table` at `scripts/tortuosity.py:418-432`,
showing that bout segmentation operates on the per-particle velocity
trace, not on FNG events:

```python
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
```

No `df_fng`, `_detect_fng_series`, or `fng_*` reference appears
anywhere in `compute_tortuosity_table` or `_segment_climbing_bouts`.

### Tests added or modified

| Test | File | Asserts |
|---|---|---|
| `test_bout_segmentation_isolates_climbing` (new) | `tests/test_tortuosity.py` | 3 climbs separated by stationary periods produce exactly 3 bouts; `bout_idx == [1, 2, 3]`. |
| `test_bout_segmentation_no_fng_dependency` (new) | `tests/test_tortuosity.py` | `compute_tortuosity_table` called with NO `df_fng` argument still produces bouts; no `event_idx` column appears. |
| `test_short_bouts_filtered` (new) | `tests/test_tortuosity.py` | A 4-frame bout below `bout_min_frames=10` is dropped, leaving only the 20-frame bout. |
| `test_low_displacement_bouts_filtered` (new) | `tests/test_tortuosity.py` | A bout with 0.75 mm net rise is dropped under `bout_min_displacement=5.0 mm`; the 15 mm bout survives. |
| `test_compute_tortuosity_table_basic_velocity_segmentation` (new, replaces `test_compute_tortuosity_table_per_event_per_particle`) | `tests/test_tortuosity.py` | Two uniformly climbing flies under `velocity_threshold=1 mm/s` produce exactly 2 rows; straight-fly T == 1, zig-zag T > 1. |
| `test_compute_tortuosity_table_no_fng_events` | DELETED | No longer meaningful (no FNG dependency to test). Replaced by `test_bout_segmentation_no_fng_dependency`. |
| `test_compute_tortuosity_table_no_particle_column` (modified) | `tests/test_tortuosity.py` | Updated to use new signature (no `df_fng`). |
| `test_compute_tortuosity_table_missing_required_column_raises` (modified) | `tests/test_tortuosity.py` | Updated to use new signature. |
| `test_detector_compute_tortuosity_writes_csv` (modified) | `tests/test_tortuosity.py` | Detector built **without** `df_fng` attribute to confirm Fix #2's FNG-independence. |

### Confirmation

**Spec item complete.**

---

## Fix #3 — Add the five config keys

### Spec requirement (verbatim)

> The original spec required five new config keys in all three ALLOWED_KEYS allowlists. The current implementation has zero new config keys; parameters are hardcoded.
>
> Requirements:
> - Add the following keys to ALL THREE ALLOWED_KEYS sites (missing any one site will silently drop the key):
>   - scripts/detector_fng.py:99-108 (load_for_gui)
>   - scripts/detector_fng.py:144-153 (load_for_main)
>   - FreeClimber_main.py:60-69 (load_parameters)
> - New keys with defaults:
>   - tortuosity_enabled = True (takes effect only when analysis_mode='individual'; allows disabling tortuosity computation)
>   - tortuosity_smoothing_window = 5 (frames; Savitzky-Golay window for Fix #4)
>   - tortuosity_velocity_threshold = 1.0 (mm/s; minimum vertical climbing velocity for Fix #2)
>   - tortuosity_bout_min_frames = 10 (minimum bout duration for Fix #2)
>   - tortuosity_bout_min_displacement = 2.0 (mm; minimum net vertical displacement for Fix #2)
> - Use the defensive getattr(self, 'key', default) pattern as in fng_* keys (see detector_fng.py:1072 and 954-957).
> - compute_tortuosity must read these from self via getattr — do NOT hardcode them.
> - The step_5 hook must guard with both flags: `if getattr(self, 'analysis_mode', 'cohort') == 'individual' and getattr(self, 'tortuosity_enabled', True):`
> - Update example/example.cfg and example_other/ex_1/climbing_1.cfg by adding the five new keys as commented-out lines showing defaults. Do not enable tortuosity in either example.
> - The status report must list each of the three ALLOWED_KEYS sites with the line numbers added, and the line range of the new getattr calls in compute_tortuosity.

### What was implemented

The five keys were appended to **every** ALLOWED_KEYS allowlist
(missing any one would silently drop the key from that load path):

| Site | File:lines (after edit) |
|---|---|
| `load_for_gui` | `scripts/detector_fng.py:102-117` — keys at lines **113-115** |
| `load_for_main` | `scripts/detector_fng.py:157-172` — keys at lines **169-171** |
| `load_parameters` | `scripts/FreeClimber_main.py:60-75` — keys at lines **72-74** |

`detector.compute_tortuosity` reads each key via defensive `getattr`
with the spec-mandated defaults at `scripts/detector_fng.py:1209-1216`:

```
1209:        # Config (Fix #3 allowlist keys, defensive getattr with safe defaults
1210:        # so any .cfg file that omits the tortuosity_* keys keeps working).
1211:        smoothing_window      = getattr(self, 'tortuosity_smoothing_window', 5)
1212:        velocity_threshold    = getattr(self, 'tortuosity_velocity_threshold', 1.0)
1213:        bout_min_frames       = getattr(self, 'tortuosity_bout_min_frames', 10)
1214:        bout_min_displacement = getattr(self, 'tortuosity_bout_min_displacement', 2.0)
1215:        pixel_to_cm           = getattr(self, 'pixel_to_cm', 1.0)
1216:        frame_rate            = getattr(self, 'frame_rate', 1.0)
```

`tortuosity_enabled` is consumed by the `step_5` guard at
`scripts/detector_fng.py:1552-1555`:

```
1552:        #---- Per-fly tortuosity metrics (individual mode + tortuosity_enabled) ----
1553:        if (getattr(self, 'analysis_mode', 'cohort') == 'individual'
1554:                and getattr(self, 'tortuosity_enabled', True)):
1555:            self.compute_tortuosity()
```

The five keys were appended as commented-out defaults to both example
configs (tortuosity NOT enabled in either):

- `example/example.cfg:44-54`
- `example_other/ex_1/climbing_1.cfg:44-54`

### Tests added or modified

No new spec-mandated tests for this fix. Coverage is exercised by the
existing detector-hook tests in `tests/test_tortuosity.py` (which
configure these keys on the synthetic detector namespace and verify
they flow through to `compute_tortuosity_table`), and by
`test_cohort_mode_unchanged` in `tests/test_linking.py` (which
verifies the `tortuosity_enabled` gate keeps cohort mode unchanged).

### Confirmation

**Spec item complete.**

---

## Fix #4 — Apply Savitzky-Golay smoothing

### Spec requirement (verbatim)

> The original spec required smoothing of (x, y) coordinates before computing path lengths to prevent noise inflation. The current implementation may not apply smoothing; there is no test_smoothing_reduces_noise_inflation in the test suite.
>
> Requirements:
> - In compute_tortuosity, apply scipy.signal.savgol_filter to x and y PER PARTICLE before computing velocities, path lengths, or bout segmentation.
> - Window from tortuosity_smoothing_window config (Fix #3); polyorder = 2.
> - Pad ends with edge values (mode='nearest' or equivalent).
> - If a particle has fewer frames than the smoothing window, skip smoothing for that particle and log a warning to stdout.
> - ALL downstream computations (velocity, path_length, net_displacement, vertical_displacement, vertical_efficiency, tortuosity, straightness, mean absolute turning angle) must use smoothed coordinates, not raw.
> - Add tests:
>   - test_smoothing_reduces_noise_inflation: a perfectly straight vertical trajectory with Gaussian noise (sigma = 0.5 mm) added to x and y produces a computed tortuosity closer to 1.0 with smoothing enabled (window=5) than without (window=1, which produces no smoothing). Assert the difference, not absolute values.
>   - test_smoothing_skipped_for_short_trajectories: a particle with 3 frames and a smoothing window of 5 produces a clear log message and computation proceeds using unsmoothed coordinates.
> - The status report must include the line range in compute_tortuosity where savgol_filter is invoked.

### What was implemented

- `from scipy.signal import savgol_filter` at
  `scripts/tortuosity.py:61`.
- Helper `_normalize_savgol_window(window)` at
  `scripts/tortuosity.py:244-262` coerces any user-supplied window to a
  positive odd integer; returns `0` to signal "skip smoothing" when
  `window <= 1`.
- Helper `smooth_xy(xs, ys, window, particle_label=None)` at
  `scripts/tortuosity.py:265-293`. **Calls
  `savgol_filter(..., polyorder=2, mode='nearest')` at lines 291-292.**
  If `n_frames < window` the function emits the stdout warning
  `!! tortuosity: skipping Savitzky-Golay smoothing for particle <id>
  (have N frames, need W); using raw coordinates` and returns the raw
  arrays unchanged.
- `compute_tortuosity_table` invokes `smooth_xy` BEFORE bout
  segmentation at `scripts/tortuosity.py:421-422`. All downstream
  metrics (`_segment_climbing_bouts` velocity trace, `path_length`,
  `net_displacement`, `vertical_displacement`, `vertical_efficiency`,
  `tortuosity`, `straightness`, `mean_absolute_turning_angle`) consume
  the smoothed arrays returned from `smooth_xy`.
- `detector.compute_tortuosity` reads `tortuosity_smoothing_window` via
  `getattr` (Fix #3) and forwards it via the `smoothing_window` kwarg
  on the call to `compute_tortuosity_table`
  (`scripts/detector_fng.py:1219-1227`).

### Tests added

| Test | File | Asserts |
|---|---|---|
| `test_smoothing_reduces_noise_inflation` | `tests/test_tortuosity.py` | Straight vertical climb + σ=0.5 mm Gaussian noise: `\|T_smoothed - 1\|` is strictly less than `\|T_raw - 1\|` (asserts the difference, not absolute values, per spec). |
| `test_smoothing_skipped_for_short_trajectories` | `tests/test_tortuosity.py` | 3-frame particle + `window=5`: captured stdout contains `"skipping Savitzky-Golay smoothing"` and `"particle 42"`; the computation still emits one bout from the raw coords. |

### Confirmation

**Spec item complete.**

---

## Full test results

Ran `python -m pytest tests/ --tb=no -q` on the final commit (c5b51fb):

| Test file | Passed | Failed |
|---|---:|---:|
| `tests/test_fng_bounds_and_detection.py` | 7 | 0 |
| `tests/test_linking.py` | 6 | 0 |
| `tests/test_tortuosity.py` | 20 | 0 |
| **Total** | **33** | **0** |

No failing tests. The 45 warnings reported are pre-existing TrackPy
warnings ("Could not generate velocity field for prediction: no
tracks") on synthetic trajectories with too few detections to fit a
predictor; unrelated to this work.

## Cohort-mode regression check

**Result: byte-identical.**

Procedure:

1. On the final commit (`c5b51fb`), ran
   `python scripts/FreeClimber_main.py --config_file example/example.cfg --process_all`
   and snapshotted `example/w1118_m_2_1.raw.csv`,
   `example/w1118_m_2_1.filtered.csv`, and
   `example/w1118_m_2_1.fng.csv` to `/tmp/postfix_*`.
2. `git stash -u && git checkout df79d0e` (the pre-fix branch tip).
3. Re-ran the same pipeline command.
4. `diff -q` between the pre-fix and post-fix outputs on all three
   files: **no differences**.
5. `git checkout feature/tortuosity-metrics && git stash pop` to
   restore working state.

Additionally, the automated regression test
`tests/test_linking.py::test_cohort_mode_unchanged` (which compares
cohort-mode `fng.csv` output byte-for-byte against a committed fixture
under both `analysis_mode='cohort'` and `analysis_mode` omitted) passes
on the post-fix branch.

Pipeline stdout in cohort mode contains no tortuosity-related lines
(no `[ Tortuosity ]` log, no `tortuosity_bouts.csv` or
`tortuosity_particle.csv` written) — the `step_5` guard
(`analysis_mode == 'individual' AND tortuosity_enabled`) short-circuits
correctly.

## Items where I considered deviating from the spec, and why I did not

1. **Particle CSV schema.** The spec explicitly names only
   `median_vertical_efficiency` as a column on
   `*.tortuosity_particle.csv`. A natural extension would be to also
   emit `median_tortuosity`, `median_straightness`, and
   `median_mean_turning_angle_rad`, since the file is otherwise sparse.
   I considered adding them but chose strict spec compliance per
   principle B (do not substitute alternatives). Per spec, the file
   carries `vial, particle, n_bouts, median_vertical_efficiency` and
   nothing else. The `n_bouts` field is included because the file
   would otherwise be ambiguous about whether `0` median means "no
   bouts" or "median of zeros".
2. **Auxiliary turning-angle metric retention.** Per principle F, the
   `mean_absolute_turning_angle` implementation and its tests were
   retained as auxiliary metrics rather than removed. They are
   explicitly described in the module docstring as auxiliary, not
   primary. No spec language requires removal.
3. **Even smoothing windows.** `scipy.signal.savgol_filter` requires
   odd window lengths. The spec specifies window=5 (odd) as the
   default. I added `_normalize_savgol_window` to round even values up
   to the next odd integer rather than raising, so a user setting
   `tortuosity_smoothing_window=6` in a cfg file gets a sensible
   `window=7` instead of a crash. This is forgiving behaviour, not a
   substitution — the spec is silent on even values. Window `<= 1`
   means "no smoothing" (which is how the
   `test_smoothing_reduces_noise_inflation` test invokes the
   no-smoothing baseline, per the spec wording "window=1, which
   produces no smoothing").
4. **`bout_idx` instead of `event_idx`.** Once bouts no longer come
   from FNG events, calling the column `event_idx` is misleading. I
   renamed to `bout_idx` (per-particle ordinal). The spec does not name
   the column.
5. **Unit-conversion explicitness.** The spec calls for "(x, y) in mm".
   I made the conversion explicit in `compute_tortuosity_table`
   (`px_to_mm = 10.0 / pixel_to_cm`) and renamed the per-bout columns
   to `path_length_mm`, `net_displacement_mm`,
   `vertical_displacement_mm` so units are visible in the output
   schema. This is documentation, not a substitution.

## Manual review items the user should check before merging

1. **Default `tortuosity_velocity_threshold = 1.0` mm/s.** This is the
   spec default. On the synthetic FNG validation fixture (where
   `pixel_to_cm=39.4`, `frame_rate=25`), a typical climbing step is
   ~50–250 mm/s, so 1.0 mm/s is well below; on real data the threshold
   may need to be tuned per rig. Recommend a quick sanity-check run on
   a representative real video before relying on the default.
2. **Default `tortuosity_bout_min_displacement = 2.0` mm.** Same caveat
   as above. With clean linked tracks this drops only very-low-rise
   bouts; on noisy data it may drop more than expected.
3. **Particle CSV columns.** As noted under "considered deviating",
   only `median_vertical_efficiency` is required. If you want median
   of tortuosity/straightness/turning angle too, that's a one-line
   extension to `compute_particle_table`.
4. **`example/w1118_m_2_1.*` outputs.** Re-running the example
   pipeline regenerated these files in place. They differ from the
   versions previously committed to git (which were stale snapshots
   from an older code version). The new files are identical between
   pre-fix and post-fix runs — i.e., my fixes don't change them — but
   they DO differ from the on-disk committed copies. Decide before
   merging whether to commit the regenerated example outputs as a
   one-time refresh. I did not commit them (work-tree only).
5. **`scripts/FreeClimber_gui.py` was NOT touched** per the spec
   exclusion. The "Per-fly tracking" GUI checkbox exists from the
   initial implementation (set `analysis_mode='individual'`) but there
   is no GUI checkbox for `tortuosity_enabled`. Anyone running from
   the GUI in individual mode will get tortuosity output with default
   parameters, which is the spec-mandated default behaviour. If a GUI
   toggle is desired later, that's a separate ticket.
6. **Smoothing edge-handling.** `mode='nearest'` is one of several
   spec-acceptable edge-padding modes ("edge values (mode='nearest' or
   equivalent)"). I chose `'nearest'` because it preserves the first
   and last point exactly, which is what `vertical_efficiency` reads
   off the bout endpoints.
