"""
Tests for per-fly trajectory linking (analysis_mode == 'individual').

Covers:
  * test_cohort_mode_unchanged           -- backward compatibility: cohort-mode
    output is byte-identical to the pre-linking behavior.
  * test_individual_mode_produces_tracks_csv -- individual mode writes a
    *.tracks.csv with a 'particle' column and no sub-threshold tracks.
  * test_per_vial_linking_isolates_vials -- per-vial linking keeps particle IDs
    disjoint across vials (no cross-vial swaps).
  * test_predictor_modes_run             -- both 'nearest_velocity' and 'none'
    predictors run without error and produce valid particle assignments.

Each test builds a minimal detector namespace from CSV/synthetic data, so no
ffmpeg or video decoding is required.
"""

import os
import sys
import types

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
FILTERED_CSV = os.path.join(_FIXTURES, "freeclimber_fng_validation_video.filtered.csv")
FNG_CSV = os.path.join(_FIXTURES, "freeclimber_fng_validation_video.fng.csv")

N_FRAMES = 750

skip_no_data = pytest.mark.skipif(
    not os.path.exists(FILTERED_CSV),
    reason=f"Validation filtered CSV not found at {FILTERED_CSV}",
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _bind(det, *names):
    """Bind the named detector methods onto a SimpleNamespace instance."""
    # Always bind the output relabeling helpers: every CSV-writing method now
    # passes its frame through _relabel_vial_col (which calls _vial_label).
    for name in tuple(names) + ("_vial_label", "_relabel_vial_col"):
        setattr(det, name, types.MethodType(getattr(dfng.detector, name), det))
    return det


def _make_det(df, tmp_dir, name="output", **cfg):
    """Build a minimal detector namespace wrapping the supplied DataFrame."""
    det = types.SimpleNamespace(
        debug=False,
        n_frames=N_FRAMES,
        frame_rate=25.0,
        pixel_to_cm=39.4,
        vials=int(df["vial"].max()) if "vial" in df.columns else 1,
        fng_enabled=True,
        fng_smooth_window=5,
        fng_climb_thresh=0.10,
        fng_fall_thresh=0.10,
        fng_min_gap=5,
        fng_recovery_thresh=0.05,
        df_filtered=df.copy(),
        file_details={},
        name_nosuffix=os.path.join(tmp_dir, name),
    )
    for key, val in cfg.items():
        setattr(det, key, val)
    return det


# ---------------------------------------------------------------------------
# test_cohort_mode_unchanged
# ---------------------------------------------------------------------------
@skip_no_data
@pytest.mark.parametrize("mode", ["cohort", None])
def test_cohort_mode_unchanged(tmp_path, mode):
    """With analysis_mode='cohort' (or omitted), output is byte-identical to the
    pre-linking behavior on the synthetic validation fixture."""
    df = pd.read_csv(FILTERED_CSV)
    cfg = {} if mode is None else {"analysis_mode": mode}
    det = _make_det(df, str(tmp_path),
                    name="freeclimber_fng_validation_video", **cfg)
    _bind(det, "_height_traces", "_detect_fng_series", "compute_fng",
          "link_trajectories")

    # Reproduce the step_5 gating: linking runs only in individual mode.
    if getattr(det, "analysis_mode", "cohort") == "individual":
        det.link_trajectories()
    det.compute_fng()

    # No linking happened -> df_filtered must not have gained a 'particle' column.
    assert "particle" not in det.df_filtered.columns

    generated = open(
        os.path.join(str(tmp_path),
                     "freeclimber_fng_validation_video.fng.csv"), "rb").read()
    reference = open(FNG_CSV, "rb").read()
    assert generated == reference, (
        "cohort-mode fng.csv is not byte-identical to the committed fixture"
    )


# ---------------------------------------------------------------------------
# test_individual_mode_produces_tracks_csv
# ---------------------------------------------------------------------------
@skip_no_data
def test_individual_mode_produces_tracks_csv(tmp_path):
    """Individual mode writes a *.tracks.csv with a 'particle' column, >= 1
    particle, and every track at least link_min_track_length frames long."""
    df = pd.read_csv(FILTERED_CSV)
    naming_keys = [c for c in df.columns
                   if c not in ("y", "x", "frame", "t", "vial")]
    file_details = {k: df[k].iloc[0] for k in naming_keys}

    det = _make_det(df, str(tmp_path),
                    name="freeclimber_fng_validation_video",
                    analysis_mode="individual",
                    link_search_range=15, link_memory=3,
                    link_predictor="nearest_velocity",
                    link_min_track_length=5)
    det.file_details = file_details
    _bind(det, "link_trajectories")
    det.link_trajectories()

    tracks_path = os.path.join(
        str(tmp_path), "freeclimber_fng_validation_video.tracks.csv")
    assert os.path.exists(tracks_path), "tracks.csv was not written"

    tracks = pd.read_csv(tracks_path)
    assert "particle" in tracks.columns, "tracks.csv missing 'particle' column"
    assert tracks["particle"].nunique() >= 1, "no particles in tracks.csv"

    lengths = tracks.groupby("particle").size()
    assert (lengths >= 5).all(), (
        f"track(s) shorter than link_min_track_length: {lengths.to_dict()}"
    )

    # Naming-convention fields propagate into tracks.csv.
    for key in naming_keys:
        assert key in tracks.columns, f"tracks.csv missing naming field '{key}'"


# ---------------------------------------------------------------------------
# test_per_vial_linking_isolates_vials
# ---------------------------------------------------------------------------
def test_per_vial_linking_isolates_vials(tmp_path):
    """Two vials, each with a fly at the SAME x, must receive disjoint particle
    IDs -- per-vial linking prevents cross-vial ID swaps."""
    frames = np.arange(30)
    rows = []
    for vial in (1, 2):
        for f in frames:
            rows.append({
                "frame": int(f),
                "t": round(f / 25.0, 3),
                "vial": vial,
                "x": 50.0,
                "y": 10.0 + float(f),
            })
    df = pd.DataFrame(rows)

    det = _make_det(df, str(tmp_path), name="two_vials",
                    analysis_mode="individual",
                    link_search_range=15, link_memory=3,
                    link_predictor="nearest_velocity",
                    link_min_track_length=5)
    _bind(det, "link_trajectories")
    det.link_trajectories()

    linked = det.df_filtered
    assert "particle" in linked.columns
    ids_v1 = set(linked.loc[linked.vial == 1, "particle"])
    ids_v2 = set(linked.loc[linked.vial == 2, "particle"])
    assert ids_v1 and ids_v2, "expected linked particles in both vials"
    assert ids_v1.isdisjoint(ids_v2), (
        f"cross-vial particle ID collision: vial 1={ids_v1} vial 2={ids_v2}"
    )


# ---------------------------------------------------------------------------
# test_predictor_modes_run
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("predictor", ["nearest_velocity", "none"])
def test_predictor_modes_run(tmp_path, predictor):
    """Both predictor modes run without error on a hand-crafted fixture of two
    flies crossing paths, and produce valid particle assignments."""
    # Fly A moves left->right, fly B moves right->left; their x-paths cross at
    # the midpoint. A small constant y-offset keeps them physically distinct.
    n = 30
    rows = []
    for f in range(n):
        t = round(f / 25.0, 3)
        rows.append({"frame": f, "t": t, "vial": 1,
                     "x": 10.0 + 2.0 * f, "y": 100.0})   # fly A
        rows.append({"frame": f, "t": t, "vial": 1,
                     "x": 70.0 - 2.0 * f, "y": 104.0})   # fly B
    df = pd.DataFrame(rows)

    det = _make_det(df, str(tmp_path), name="crossing_" + predictor,
                    analysis_mode="individual",
                    link_search_range=15, link_memory=3,
                    link_predictor=predictor,
                    link_min_track_length=5)
    _bind(det, "link_trajectories")
    det.link_trajectories()

    linked = det.df_filtered
    assert not linked.empty, f"predictor '{predictor}' produced no tracks"
    assert "particle" in linked.columns

    lengths = linked.groupby("particle").size()
    # Two flies present in every frame -> two full-length tracks expected.
    assert len(lengths) >= 2, (
        f"predictor '{predictor}': expected >=2 tracks, got {len(lengths)}"
    )
    assert (lengths >= 5).all(), (
        f"predictor '{predictor}': stub track survived filter_stubs: "
        f"{lengths.to_dict()}"
    )
