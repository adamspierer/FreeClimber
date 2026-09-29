"""
Regression tests for audit bug fixes in detector_fng.py and FreeClimber_main.py.

Covers:
  * config parsing keeps Windows paths intact (no '\\t' -> tab mangling)
  * file names with 3-letter extensions (.mov/.mp4) are not truncated
  * window is coerced to an int when the crop range is shorter than it
  * invert_y measures height from the ROI floor, not the lowest detection
  * local_linear_regression / get_slopes tolerate an empty vial
  * image_plot draws a frame that has no detections (used to crash step 6)
  * vials.txt 'n =' fly-count line and bare comma-separated ID lists
  * FreeClimber_main.file_walker matches the file suffix case-insensitively

No ffmpeg or video decoding is required.
"""

import os
import sys
import types

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "scripts"))
import detector_fng as dfng  # noqa: E402
import FreeClimber_main as fcm  # noqa: E402


def _bind(det, *names):
    for name in names:
        setattr(det, name, types.MethodType(getattr(dfng.detector, name), det))
    return det


def test_parse_config_line_windows_path():
    key, val = dfng.parse_config_line(r'path_project="D:\test\new\run"')
    assert (key, val) == ("path_project", r"D:\test\new\run")
    assert dfng.parse_config_line("background_image=None") == ("background_image", None)
    assert dfng.parse_config_line("vials=5") == ("vials", 5)
    assert dfng.parse_config_line("flies_per_vial=[10, 9]") == ("flies_per_vial", [10, 9])
    assert dfng.parse_config_line("not_a_key=1") is None
    assert dfng.parse_config_line("# comment") is None


def test_specify_paths_three_letter_extension():
    det = types.SimpleNamespace(debug=False, path_project="x",
                                naming_convention="geno_sex_day_rep",
                                vial_id_vars=4, vials=2, vial_color_map=plt.cm.jet)
    _bind(det, "specify_paths_details")
    det.specify_paths_details(os.path.join("folder", "w1118_m_2_1.mov"))
    assert det.name == "w1118_m_2_1"
    assert det.name_nosuffix == os.path.join("folder", "w1118_m_2_1")
    assert det.file_details == {"geno": "w1118", "sex": "m", "day": "2", "rep": "1"}
    assert det.vial_ID == ["w1118", "m", "2", "1"]  # last field no longer blanked


def test_window_becomes_int():
    det = types.SimpleNamespace(debug=False, vials=1, diameter=7, frame_rate=30,
                                blank_0=0, blank_n=40, crop_0=0, crop_n=40,
                                window=50, check_frame=100)
    _bind(det, "check_variable_formats")
    det.check_variable_formats()
    assert det.window == 32 and isinstance(det.window, int)
    assert det.check_frame == 39  # crop_n is exclusive


def test_invert_y_uses_roi_floor():
    det = types.SimpleNamespace(debug=False, h=300)
    _bind(det, "invert_y", "_floor_px")
    spots = pd.DataFrame({"y": [290.0, 100.0]})
    assert det.invert_y(spots).tolist() == [10.0, 200.0]
    det.floor_y = 280
    assert det.invert_y(spots).tolist() == [-10.0, 180.0]


def _slopes_det(df):
    det = types.SimpleNamespace(debug=False, vials=2, vial_labels=None, crop_0=0,
                                crop_n=60, window=20, vial_ID=["a", "b"], df_filtered=df)
    return _bind(det, "local_linear_regression", "get_slopes", "_vial_label")


def test_get_slopes_tolerates_empty_vial():
    frames = np.arange(60)
    df = pd.DataFrame({"frame": frames, "vial": 2, "y": 3.0 * frames})
    det = _slopes_det(df)
    assert det.local_linear_regression(df.iloc[0:0]).empty

    det.get_slopes()
    assert np.isnan(det.result[1][1])        # empty vial 1: NaN row, not skipped
    assert det.result[1][0] == "a_b_1"
    assert det.result[2][3] == 3.0           # vial 2 slope
    assert det.result[3][0] == "a_b_all"


def test_image_plot_frame_without_spots():
    det = types.SimpleNamespace(debug=False, clean_stack=np.zeros((10, 20, 30)),
                                h=20, w=30, bin_lines=[0, 15, 30], vials=2,
                                vial_color_map=plt.cm.jet)
    _bind(det, "image_plot")
    spots = pd.DataFrame({"frame": [1, 2], "x": [5.0, 6.0], "y": [5.0, 6.0], "vial": [1, 2]})
    fig, ax = plt.subplots()
    det.image_plot(df=spots, frame=7, ax=ax)   # frame 7 has no spots
    assert "no spots" in ax.get_title()
    det.image_plot(df=spots, frame=10, ax=ax)  # past the end -> last frame
    plt.close(fig)


def test_vials_sidecar_fly_counts_and_bare_ids(tmp_path):
    det = types.SimpleNamespace(debug=False, vials=5, vial_labels=None)
    _bind(det, "apply_vials_sidecar")
    det._parse_vial_labels = dfng.detector._parse_vial_labels
    (tmp_path / "vials.txt").write_text("b26, b27, b29\nn = 10, 9, 10\n")
    det.apply_vials_sidecar(str(tmp_path / "clip.mov"))
    assert det.vial_labels == ["b26", "b27", "b29"] and det.vials == 3
    assert det.flies_per_vial == [10, 9, 10]

    (tmp_path / "vials.txt").write_text("flies: 8\n")
    det.apply_vials_sidecar(str(tmp_path / "clip.mov"))
    assert det.flies_per_vial == 8


def test_file_walker_case_insensitive(tmp_path):
    for name in ("a.MOV", "b.mov", "c.slopes.csv", "b.slopes.csv"):
        (tmp_path / name).write_text("")
    fc = object.__new__(fcm.FreeClimber)
    fc.args = types.SimpleNamespace(debug=False)
    found = fc.file_walker(str(tmp_path), endswith="mov")
    assert [os.path.basename(f) for f in found] == ["a.MOV", "b.mov"]
    undone = fc.file_walker(str(tmp_path), endswith="mov", undone=True)
    assert [os.path.basename(f) for f in undone] == ["a.MOV"]
