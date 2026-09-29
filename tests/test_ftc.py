"""
Tests for the failure-to-climb (FTC) measure (detector.compute_ftc) and the
FNG noise floor that keeps non-climbing vials from being scored as falls.

Covers:
  * test_ftc_cohort_detected_counts   -- cohort mode, no fly count configured:
    FTC is the flies below the line out of those detected.
  * test_ftc_expected_counts_and_warning -- with flies_per_vial, FTC counts
    flies that never crossed the line even when they are not detected, and a
    count_warning flags the shortfall.
  * test_ftc_particle_outcomes        -- individual mode: climber / fng / ftc /
    unscored outcomes and latency per fly.
  * test_ftc_window                   -- ftc_start_frame / ftc_window_sec limit
    the assessment window.
  * test_ftc_disabled                 -- ftc_enabled=False writes nothing.
  * test_fng_ignores_non_climbing_vial -- a vial of motionless, jittery flies
    produces no FNG events (it used to produce many).

Synthetic data only; no ffmpeg or video decoding is required.
"""

import os
import sys
import types

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "scripts"))
import detector_fng as dfng  # noqa: E402

PX_PER_CM = 20.0
FPS = 10.0
N_FRAMES = 100


def _bind(det):
    for name in ("compute_ftc", "_ftc_window", "_flies_expected", "_floor_px",
                 "_ftc_line_px", "_detect_fng_series", "_vial_label",
                 "_relabel_vial_col"):
        setattr(det, name, types.MethodType(getattr(dfng.detector, name), det))
    return det


def _make_det(df, tmp_path, **cfg):
    det = types.SimpleNamespace(
        debug=False, vials=int(df.vial.max()), vial_labels=None,
        pixel_to_cm=PX_PER_CM, frame_rate=FPS, n_frames=N_FRAMES,
        crop_0=0, crop_n=N_FRAMES,
        fng_smooth_window=5, fng_climb_thresh=0.10, fng_fall_thresh=0.10,
        fng_min_gap=5,
        ftc_height_cm=2.0, ftc_eval_frames=5, ftc_min_coverage=0.8,
        df_filtered=df, file_details={"geno": "w1118"},
        name_nosuffix=str(tmp_path / "clip"),
    )
    for key, val in cfg.items():
        setattr(det, key, val)
    return _bind(det)


def _track(vial, particle, heights_cm, frames=None, jitter=0.0, seed=0):
    frames = np.arange(len(heights_cm)) if frames is None else np.asarray(frames)
    rng = np.random.default_rng(seed)
    y = np.asarray(heights_cm, dtype=float) * PX_PER_CM + rng.normal(0, jitter, len(frames))
    return pd.DataFrame({"frame": frames, "vial": vial, "particle": particle,
                         "x": 10.0 * vial, "y": y})


def _climb(to_cm, over=30, n=N_FRAMES):
    """Rise linearly from 0.1 cm to to_cm over 'over' frames, then hold."""
    return np.concatenate([np.linspace(0.1, to_cm, over), np.full(n - over, to_cm)])


def _scenario():
    """vial 1: 3 climbers; vial 2: 3 motionless flies; vial 3: one fly that
    climbs above the line then falls, one partial climber (1 cm)."""
    parts = []
    for p in range(3):
        parts.append(_track(1, 10 + p, _climb(5.0), jitter=0.3, seed=p))
    for p in range(3):
        parts.append(_track(2, 20 + p, np.full(N_FRAMES, 0.2), jitter=0.5, seed=10 + p))
    up_then_down = np.concatenate([np.linspace(0.1, 5.0, 30), np.full(20, 5.0),
                                   np.linspace(5.0, 0.3, 5), np.full(45, 0.3)])
    parts.append(_track(3, 30, up_then_down, jitter=0.2, seed=20))
    parts.append(_track(3, 31, _climb(1.0), jitter=0.2, seed=21))
    return pd.concat(parts, ignore_index=True)


def test_ftc_cohort_detected_counts(tmp_path):
    df = _scenario().drop(columns="particle")
    det = _make_det(df, tmp_path)
    det.compute_ftc()

    out = det.df_ftc.set_index("vial")
    assert list(out.method.unique()) == ["detected"]
    assert out.loc[1, "ftc_count"] == 0 and out.loc[1, "n_reached_line"] == 3
    assert out.loc[2, "ftc_count"] == 3 and out.loc[2, "ftc_fraction"] == 1.0
    # the faller reached the line before falling: only the partial climber failed
    assert out.loc[3, "ftc_count"] == 1 and out.loc[3, "n_above_line_end"] == 0
    assert (out.count_warning == "").all()
    assert out["n_tracks_ftc"].isna().all()  # no per-fly counts in cohort mode

    written = pd.read_csv(str(tmp_path / "clip.ftc.csv"))
    assert written.columns[0] == "geno"  # naming-convention details prepended
    assert not os.path.exists(str(tmp_path / "clip.ftc_particle.csv"))


def test_ftc_expected_counts_and_warning(tmp_path):
    df = _scenario().drop(columns="particle")
    # 5 flies loaded per vial, but only 3 (or 2) are ever detected
    det = _make_det(df, tmp_path, flies_per_vial=[3, 5, 2])
    det.compute_ftc()

    out = det.df_ftc.set_index("vial")
    assert list(out.method.unique()) == ["expected"]
    assert out.loc[1, "ftc_count"] == 0 and out.loc[1, "count_warning"] == ""
    # the two undetected flies in vial 2 count as failures
    assert out.loc[2, "ftc_count"] == 5 and out.loc[2, "ftc_fraction"] == 1.0
    assert out.loc[2, "count_warning"] == "fewer_detected_than_expected"


def test_ftc_climbers_that_vanish_still_counted(tmp_path):
    """Climbers that stop at the top and drop out of the detections (subtracted
    into the background) still count as having reached the line, and the
    vial is flagged."""
    df = _scenario().drop(columns="particle")
    df = df[~((df.vial == 1) & (df.frame >= 50))]  # vial 1 flies vanish at frame 50
    det = _make_det(df, tmp_path)
    det.compute_ftc()
    row = det.df_ftc.set_index("vial").loc[1]
    assert row.n_reached_line == 3 and row.ftc_count == 0
    assert row.n_detected_end == 0 and row.count_warning == "detections_dropped"


def test_ftc_particle_outcomes(tmp_path):
    df = _scenario()
    partial_fall = np.concatenate([np.linspace(0.1, 1.2, 25), np.full(10, 1.2),
                                   np.linspace(1.2, 0.2, 5), np.full(60, 0.2)])
    df = pd.concat([df,
                    _track(3, 32, np.full(10, 0.2), frames=np.arange(10)),  # brief, low -> unscored
                    _track(3, 33, partial_fall, jitter=0.2, seed=22)],     # partial climb, fell -> fng
                   ignore_index=True)
    det = _make_det(df, tmp_path)
    det.compute_ftc()

    parts = det.df_ftc_particle.set_index("particle")
    assert (parts.loc[[10, 11, 12], "outcome"] == "climber").all()
    assert (parts.loc[[20, 21, 22], "outcome"] == "ftc").all()
    assert parts.loc[30, "outcome"] == "fng"
    assert parts.loc[31, "outcome"] == "ftc"  # partial climb never reached 2 cm
    assert parts.loc[32, "outcome"] == "unscored"
    # a partial climb that ends in a fall is a fall, even below the line
    assert parts.loc[33, "outcome"] == "fng"
    assert not parts.loc[33, "reached_line"] and parts.loc[33, "n_falls"] >= 1

    # climbers reach 2 cm about 11-12 frames into a 30-frame climb to 5 cm
    assert parts.loc[10, "latency_sec"] == pytest.approx(1.15, abs=0.2)
    assert parts.loc[[20, 31], "latency_sec"].isna().all()

    counts = det.df_ftc.set_index("vial")
    assert counts.loc[1, "n_tracks_climber"] == 3
    assert counts.loc[2, "n_tracks_ftc"] == 3
    assert counts.loc[3, ["n_tracks_fng", "n_tracks_ftc", "n_tracks_unscored"]].tolist() == [2, 1, 1]
    assert os.path.exists(str(tmp_path / "clip.ftc_particle.csv"))


def test_ftc_window(tmp_path):
    df = _scenario()
    # 1 s window starting at frame 0: nobody has reached 2 cm yet
    det = _make_det(df, tmp_path, ftc_window_sec=1.0)
    det.compute_ftc()
    row = det.df_ftc.iloc[0]
    assert (row.window_start_frame, row.window_end_frame, row.window_sec) == (0, 9, 1.0)
    assert (det.df_ftc_particle.reached_line == False).all()  # noqa: E712

    # ftc_start_frame is an absolute video frame; frames are crop-relative
    det = _make_det(df, tmp_path, crop_0=20, crop_n=120, ftc_start_frame=50)
    assert det._ftc_window(df) == (30, 99)


def test_ftc_line_defaults_to_roi_top(tmp_path):
    """With no ftc_height_cm the line is the top of the ROI box, one spot
    diameter down (TrackPy cannot detect spots centred on the edge)."""
    df = _scenario()
    det = _make_det(df, tmp_path, ftc_height_cm=None, h=110, diameter=7)
    assert det._ftc_line_px() == 103.0
    det.compute_ftc()
    assert det.df_ftc.ftc_height_cm.iloc[0] == pytest.approx(103.0 / PX_PER_CM)
    # climbers top out at 5 cm (100 px): below a 5.15 cm line, so none reach it
    assert not det.df_ftc_particle.reached_line.any()

    det = _make_det(df, tmp_path, ftc_height_cm=None, h=110, floor_y=100, diameter=7)
    assert det._ftc_line_px() == 93.0  # measured from the configured floor


def test_ftc_disabled(tmp_path):
    det = _make_det(_scenario(), tmp_path, ftc_enabled=False)
    det.compute_ftc()
    assert det.df_ftc.empty
    assert not os.path.exists(str(tmp_path / "clip.ftc.csv"))


def test_fng_ignores_non_climbing_vial():
    """Motionless flies (sub-millimetre jitter) must not be scored as falls:
    the 0..1 normalisation used to stretch the jitter into fake climb/falls."""
    rng = np.random.default_rng(0)
    jitter = pd.Series(5 + rng.normal(0, 1.5, 300), index=range(300))
    det = types.SimpleNamespace(pixel_to_cm=26.0, frame_rate=30.0, n_frames=300)
    detect = types.MethodType(dfng.detector._detect_fng_series, det)

    assert len(detect(jitter)) == 0
    det.fng_min_range_cm = 0.0  # the old behaviour, for comparison
    assert len(detect(jitter)) > 5
