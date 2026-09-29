FreeClimber-FNG Test Sample Video
==================================

PURPOSE
  Synthetic single-vial, single-fly video designed to test detection
  parameters and FNG event-finding logic against known ground truth.

VIDEO PROPERTIES
  resolution      : 480 x 720 px
  frame rate      : 25 fps
  duration        : 30 s (750 frames)
  background gray : 235
  fly intensity   : Gaussian, peak=20, sigma=2.0 px
                    (~9 px FWHM diameter)

VIAL ROI (pixel coordinates, image-space)
  x range : 200 to 280
  y range : 80 (top) to 660 (bottom)
  inner height : 580 px
  physical scale assumed : 9.5 cm

SUGGESTED FREECLIMBER PARAMETERS
  diameter        : 9 (odd integer, matches fly FWHM)
  threshold       : low (e.g. 1-3) — fly is high-contrast on flat bg
  minmass         : start ~50, raise if false positives appear
  maxsize         : leave default
  ecc_low / ecc_high : 0 / 1 (fly is circular by construction)
  invert          : True (dark fly on light background)
  vials           : 1
  vial_id         : 1
  diameter expected : single column near x=240

GROUND-TRUTH FALL EVENTS (5 total)
  Event 1 : peak frame 150, fall 151-165 (14 frames),
            rise_norm 0.92 → end 0.40, drop 0.52 (~5.0 cm)
  Event 2 : peak frame 280, fall 281-293 (13 frames),
            rise_norm 0.85 → end 0.55, drop 0.30 (~2.9 cm)
  Event 3 : peak frame 410, fall 411-428 (18 frames),
            rise_norm 0.95 → end 0.30, drop 0.65 (~6.2 cm)
  Event 4 : peak frame 540, fall 541-555 (15 frames),
            rise_norm 0.78 → end 0.50, drop 0.28 (~2.7 cm)
  Event 5 : peak frame 660, fall 661-680 (20 frames),
            rise_norm 0.90 → end 0.20, drop 0.70 (~6.7 cm)

NOTES ON DETECTION DIFFICULTY
  This sample is deliberately easy — high contrast, no occlusion,
  smooth trajectory with small Gaussian jitter, no spot identity
  ambiguity (one fly), no events shorter than 12 frames. Successful
  detection here confirms the pipeline runs end-to-end and that
  parameters are roughly tuned. It does NOT stress-test the algorithm
  for the harder cases that arise in real videos: multiple flies in
  one vial, brief identity swaps, low-position wandering, sub-5-frame
  events, or wall-clinging behavior.
