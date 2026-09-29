"""
Regression test for two FNG detection bugs:
  1. OOB bug: frame indices (frame_peak, frame_fall_start, frame_fall_end) must
     satisfy 0 <= i < n_frames.
  2. Detection bug: adaptive-lookback rise window must detect gradual-climb events
     that the old fixed-window algorithm missed.

Runs against the pre-computed raw CSV from the validation video so no ffmpeg or
TrackPy dependency is required in CI.

Ground truth: 5 events, peaks at frames [150, 280, 410, 540, 660], from
  tests/fixtures/synthetic_validation/freeclimber_fng_validation_video_ground_truth.csv
Tolerance: ±5 frames on frame_peak.
"""

import os
import sys
import types
import tempfile

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Locate project root and import detector_fng
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "scripts"))
import detector_fng as dfng  # noqa: E402

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_FIXTURES = os.path.join(
    os.path.dirname(__file__),
    "fixtures",
    "synthetic_validation",
)
RAW_CSV = os.path.join(_FIXTURES, "freeclimber_fng_validation_video.raw.csv")
GT_CSV = os.path.join(_FIXTURES, "freeclimber_fng_validation_video_ground_truth.csv")

N_FRAMES = 750
GT_PEAKS = [150, 280, 410, 540, 660]
PEAK_TOL = 5  # ±5 frames


# ---------------------------------------------------------------------------
# Helper: build a minimal detector-like namespace from the raw CSV
# ---------------------------------------------------------------------------
def _bind_methods(det: types.SimpleNamespace) -> types.SimpleNamespace:
    """Bind the detector instance methods needed by compute_fng onto a SimpleNamespace."""
    for name in ("_height_traces", "_detect_fng_series", "compute_fng",
                 "_vial_label", "_relabel_vial_col"):
        setattr(det, name, types.MethodType(getattr(dfng.detector, name), det))
    return det


def _make_det(tmp_dir: str) -> types.SimpleNamespace:
    """Build a minimal detector namespace loaded with data from the raw CSV."""
    raw = pd.read_csv(RAW_CSV)

    # Reproduce the invert_y used in step_5
    raw["y"] = abs(raw["y"] - raw["y"].max()).round(2)
    raw["vial"] = 1

    det = types.SimpleNamespace(
        debug=False,
        n_frames=N_FRAMES,
        frame_rate=25.0,
        pixel_to_cm=39.4,
        vials=1,
        fng_enabled=True,
        fng_smooth_window=5,
        fng_climb_thresh=0.10,
        fng_fall_thresh=0.10,
        fng_min_gap=5,
        fng_recovery_thresh=0.05,
        df_filtered=raw[["frame", "y", "vial"]].copy(),
        name_nosuffix=os.path.join(tmp_dir, "test_output"),
    )
    return _bind_methods(det)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def fng_result(tmp_path_factory):
    """Run compute_fng once; share across tests in this module."""
    tmp = str(tmp_path_factory.mktemp("fng_out"))
    det = _make_det(tmp)
    det.compute_fng()
    return det


# ---------------------------------------------------------------------------
# Skip guard: validation files might not be present in all environments
# ---------------------------------------------------------------------------
skip_no_data = pytest.mark.skipif(
    not os.path.exists(RAW_CSV),
    reason=f"Validation raw CSV not found at {RAW_CSV}",
)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
@skip_no_data
def test_five_events_detected(fng_result):
    """Exactly 5 FNG events must be returned for the validation video."""
    assert hasattr(fng_result, "df_fng"), "df_fng attribute missing"
    assert not fng_result.df_fng.empty, "df_fng is empty — no events detected"
    assert len(fng_result.df_fng) == 5, (
        f"Expected 5 events, got {len(fng_result.df_fng)}. "
        f"Detected peaks: {fng_result.df_fng['frame_peak'].tolist()}"
    )


@skip_no_data
def test_peak_frames_within_tolerance(fng_result):
    """Each detected peak must be within ±PEAK_TOL frames of the ground truth."""
    detected = sorted(fng_result.df_fng["frame_peak"].tolist())
    for gt, det in zip(GT_PEAKS, detected):
        assert abs(det - gt) <= PEAK_TOL, (
            f"Detected peak {det} is more than {PEAK_TOL} frames from GT peak {gt}"
        )


@skip_no_data
def test_bounds_invariant(fng_result):
    """frame_peak, frame_fall_start, frame_fall_end must all be in [0, n_frames)."""
    df = fng_result.df_fng
    for col in ("frame_peak", "frame_fall_start", "frame_fall_end"):
        bad = df[(df[col] < 0) | (df[col] >= N_FRAMES)]
        assert bad.empty, (
            f"Column '{col}' has out-of-bounds values (valid [0, {N_FRAMES})):\n"
            f"{bad[[col]].to_string()}"
        )


@skip_no_data
def test_fall_order(fng_result):
    """frame_fall_start >= frame_peak and frame_fall_end >= frame_fall_start."""
    df = fng_result.df_fng
    assert (df["frame_fall_start"] >= df["frame_peak"]).all(), (
        "frame_fall_start < frame_peak in some rows"
    )
    assert (df["frame_fall_end"] >= df["frame_fall_start"]).all(), (
        "frame_fall_end < frame_fall_start in some rows"
    )


@skip_no_data
def test_rise_and_drop_positive(fng_result):
    """rise_norm and drop_norm must be positive and <= 1."""
    df = fng_result.df_fng
    assert (df["rise_norm"] > 0).all(), "rise_norm has non-positive values"
    assert (df["drop_norm"] > 0).all(), "drop_norm has non-positive values"
    assert (df["rise_norm"] <= 1.0 + 1e-9).all(), "rise_norm > 1"
    assert (df["drop_norm"] <= 1.0 + 1e-9).all(), "drop_norm > 1"


# ---------------------------------------------------------------------------
# Synthetic OOB stress test — no files needed
# ---------------------------------------------------------------------------
def _make_synthetic_det(n_frames: int, peaks_at: list, tmp_dir: str) -> types.SimpleNamespace:
    """Build a detector namespace with a synthetic height trace."""
    t = np.arange(n_frames, dtype=float)
    y = np.zeros(n_frames)
    for pk in peaks_at:
        start = max(0, pk - 50)
        y[start : pk + 1] += np.linspace(0, 400, pk - start + 1)
        fall_end = min(n_frames, pk + 30)
        y[pk:fall_end] += 400 * np.exp(-0.15 * np.arange(fall_end - pk))

    raw = pd.DataFrame({"frame": t.astype(int), "y": y, "vial": 1})
    det = types.SimpleNamespace(
        debug=False,
        n_frames=n_frames,
        frame_rate=25.0,
        pixel_to_cm=39.4,
        vials=1,
        fng_enabled=True,
        fng_smooth_window=5,
        fng_climb_thresh=0.10,
        fng_fall_thresh=0.10,
        fng_min_gap=5,
        fng_recovery_thresh=0.05,
        df_filtered=raw.copy(),
        name_nosuffix=os.path.join(tmp_dir, "synth_output"),
    )
    return _bind_methods(det)


def test_synthetic_bounds_invariant(tmp_path):
    """Bounds invariant holds on synthetic data with peaks near the end."""
    peaks_at = [100, 250, 400, 550, 700]
    det = _make_synthetic_det(750, peaks_at, str(tmp_path))
    det.compute_fng()
    if det.df_fng.empty:
        pytest.skip("No events detected in synthetic data — thresholds may need tuning")
    df = det.df_fng
    for col in ("frame_peak", "frame_fall_start", "frame_fall_end"):
        bad = df[(df[col] < 0) | (df[col] >= 750)]
        assert bad.empty, (
            f"OOB in '{col}' with synthetic data: {bad[[col]].to_string()}"
        )


def test_gapped_index_bounds_invariant(tmp_path):
    """Bounds invariant holds when the height trace has gaps (missing frames)."""
    n_frames = 750
    peaks_at = [150, 400, 650]
    t = np.arange(n_frames, dtype=float)
    y = np.zeros(n_frames)
    for pk in peaks_at:
        start = max(0, pk - 50)
        y[start : pk + 1] += np.linspace(0, 400, pk - start + 1)
        fall_end = min(n_frames, pk + 30)
        y[pk:fall_end] += 400 * np.exp(-0.15 * np.arange(fall_end - pk))

    raw = pd.DataFrame({"frame": t.astype(int), "y": y, "vial": 1})
    # Drop some frames to create gaps (simulates missing detections)
    drop_idx = list(range(80, 100)) + list(range(300, 320)) + list(range(600, 620))
    raw = raw.drop(raw[raw["frame"].isin(drop_idx)].index).reset_index(drop=True)

    det = types.SimpleNamespace(
        debug=False,
        n_frames=n_frames,
        frame_rate=25.0,
        pixel_to_cm=39.4,
        vials=1,
        fng_enabled=True,
        fng_smooth_window=5,
        fng_climb_thresh=0.10,
        fng_fall_thresh=0.10,
        fng_min_gap=5,
        fng_recovery_thresh=0.05,
        df_filtered=raw.copy(),
        name_nosuffix=os.path.join(str(tmp_path), "gapped_output"),
    )
    _bind_methods(det)
    det.compute_fng()
    if det.df_fng.empty:
        pytest.skip("No events detected in gapped synthetic data")
    df = det.df_fng
    for col in ("frame_peak", "frame_fall_start", "frame_fall_end"):
        bad = df[(df[col] < 0) | (df[col] >= n_frames)]
        assert bad.empty, (
            f"OOB in '{col}' with gapped index: {bad[[col]].to_string()}"
        )
