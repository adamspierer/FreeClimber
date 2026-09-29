"""
Tests for the per-folder ``vials.txt`` sidecar (detector.apply_vials_sidecar)
and the vial-ID labeling it enables.

Covers:
  * test_count_only            -- a bare integer / ``vials=`` line overrides the
    vial count (original backward-compatible behavior) and sets no labels.
  * test_id_line_sets_count_and_labels -- an ``id =`` line sets both the labels
    and the vial count from the number of IDs listed.
  * test_id_count_mismatch_id_wins -- when a count line disagrees with the
    number of IDs, the ID count wins.
  * test_missing_and_empty_noop -- a missing or comment-only file leaves the
    cfg vials unchanged and labels unset.
  * test_vial_label_mapping     -- _vial_label maps positional indices to IDs
    and passes through out-of-range / 'all' values.
  * test_relabel_vial_col       -- _relabel_vial_col rewrites the 'vial' column
    of an output DataFrame without mutating the input or touching other columns.

Each test binds the relevant methods onto a minimal namespace, so no ffmpeg or
video decoding is required.
"""

import os
import sys
import types

import pandas as pd
import pytest

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "scripts"))
import detector_fng as dfng  # noqa: E402


def _bind(det, *names):
    for name in names:
        setattr(det, name, types.MethodType(getattr(dfng.detector, name), det))
    return det


def _make_det(vials=5):
    det = types.SimpleNamespace(debug=False, vials=vials, vial_labels=None)
    _bind(det, "apply_vials_sidecar", "_vial_label", "_relabel_vial_col")
    # _parse_vial_labels is a @staticmethod -- attach it unbound (no self).
    det._parse_vial_labels = staticmethod(dfng.detector._parse_vial_labels).__func__
    return det


def _sidecar(tmp_path, text):
    """Write vials.txt into tmp_path and return a fake video path beside it."""
    (tmp_path / "vials.txt").write_text(text)
    return str(tmp_path / "clip.mp4")


def test_count_only(tmp_path):
    det = _make_det(vials=5)
    det.apply_vials_sidecar(_sidecar(tmp_path, "# folder of 3-vial videos\n3\n"))
    assert det.vials == 3
    assert det.vial_labels is None


@pytest.mark.parametrize("text", [
    "id = 33, 34, 35\n",
    "ids: 33,34,35\n",
    "# kept the right three\nid = 33, 34, 35\n",
])
def test_id_line_sets_count_and_labels(tmp_path, text):
    det = _make_det(vials=5)
    det.apply_vials_sidecar(_sidecar(tmp_path, text))
    assert det.vial_labels == [33, 34, 35]
    assert det.vials == 3  # count derived from the number of IDs


def test_id_count_mismatch_id_wins(tmp_path):
    det = _make_det(vials=5)
    det.apply_vials_sidecar(_sidecar(tmp_path, "vials = 5\nid = 33, 34, 35\n"))
    assert det.vial_labels == [33, 34, 35]
    assert det.vials == 3


def test_non_numeric_labels_allowed(tmp_path):
    det = _make_det(vials=5)
    det.apply_vials_sidecar(_sidecar(tmp_path, "id = 33, B34, 35\n"))
    assert det.vial_labels == [33, "B34", 35]
    assert det.vials == 3


@pytest.mark.parametrize("text", [None, "", "# nothing but a comment\n"])
def test_missing_and_empty_noop(tmp_path, text):
    det = _make_det(vials=5)
    if text is None:
        video = str(tmp_path / "clip.mp4")  # no vials.txt written
    else:
        video = _sidecar(tmp_path, text)
    det.apply_vials_sidecar(video)
    assert det.vials == 5
    assert det.vial_labels is None


def test_vial_label_mapping():
    det = _make_det()
    det.vial_labels = [33, 34, 35]
    assert [det._vial_label(i) for i in (1, 2, 3)] == [33, 34, 35]
    assert det._vial_label(4) == 4        # out of range -> unchanged
    assert det._vial_label("all") == "all"  # cohort aggregate -> unchanged

    det.vial_labels = None
    assert det._vial_label(2) == 2        # no labels -> positional


def test_relabel_vial_col():
    det = _make_det()
    det.vial_labels = [33, 34, 35]
    df = pd.DataFrame({"vial": [1, 2, 3, 1], "x": [10, 20, 30, 40]})

    out = det._relabel_vial_col(df)
    assert list(out["vial"]) == [33, 34, 35, 33]
    assert list(out["x"]) == [10, 20, 30, 40]      # other columns untouched
    assert list(df["vial"]) == [1, 2, 3, 1]        # input not mutated

    # No labels configured -> returned unchanged.
    det.vial_labels = None
    assert det._relabel_vial_col(df) is df
