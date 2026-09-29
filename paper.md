---
title: "FreeClimber-FNG: Automated Detection of Failed Negative Geotaxis and Fall Distance in Drosophila Climbing Assays"
tags:
  - Python
  - Drosophila
  - behavioral analysis
  - reproducible research
authors:
  - name: Jordan Vasu
    orcid: 0009-0001-3178-4229
    affiliation: "1"
  - name: Richard E. Hartman
    orcid: 0000-0001-9235-3169
    affiliation: "1"
affiliations:
  - name: Loma Linda University, School of Behavioral Health, Department of Psychology
    index: "1"
date: 27 December 2025
bibliography: paper.bib
---

# Summary

Negative geotaxis (climbing) assays in *Drosophila melanogaster* are widely used to classify motor performance and subsequent disruption by age, traumatic brain injury, genetic perturbation, or pharmacological manipulation. The open-source *FreeClimber* platform automates extraction of climbing velocity from video recordings using background subtraction, particle tracking, and local linear regression on vial-wise mean vertical position traces [@spierer2021freeclimber].

Investigators can additionally score failed negative geotaxis (FNG) events. FNG is defined as an event in which a fly climbs but then loses footing and falls within the assay window. FNG is typically hand-scored, which is time-consuming, rater-dependent, and difficult to reproduce.

*FreeClimber-FNG* is a Python extension to *FreeClimber* that adds automated detection of FNG events and computation of fall distance using the same frame-resolved position data already produced by *FreeClimber*. The extension runs automatically during the standard pipeline, requires no GUI changes, and yields machine-readable `.fng.csv` output suitable for statistical analysis. This tool provides a reproducible, scriptable alternative to manual FNG scoring.

# Statement of Need

*Drosophila* climbing assays remain a central behavioral measure in neuroscience and aging research due to their simplicity, experimental flexibility, and sensitivity to subtle motor impairments [@gargano2005ring]. Laboratories may choose to report the presence of FNG events in addition to standard velocity and height metrics, as falls may characterize motor discoordination, intoxication, fatigue, or insult-related deficit. However:

1. No commonly used open-source tool automates FNG detection
2. Most laboratories rely on manual review of videos
3. Scoring is difficult to standardize across raters
4. Results are not easily auditable or reproducible

*FreeClimber-FNG* addresses this gap by integrating directly with *FreeClimber's* existing data pipeline and by adding an analysis layer for identifying climb–fall sequences and quantifying fall distance. Researchers who already use *FreeClimber* can therefore obtain both velocity, height, and fall-behavior measures from the same recorded video without additional steps.

# Research Impact Statement

FreeClimber-FNG is archived as a citable software deposit on Zenodo [@vasu2025freeclimberfng], providing a versioned, DOI-linked release suitable for citation in methods sections. The software is currently in active use in our laboratory, where it has been applied to climbing assays across wild-type, traumatic brain injury, and pharmacological manipulation conditions. Validation manuscripts reporting these datasets are in preparation and will cite the software directly upon posting as preprints.

FreeClimber, the upstream platform this tool extends, has been adopted across multiple laboratories and published in a peer-reviewed methods paper [@spierer2021freeclimber], providing an established ecosystem into which FreeClimber-FNG integrates.

# Software Design

FreeClimber-FNG is designed as a modular extension to the existing FreeClimber analysis pipeline. The software reuses FreeClimber's established video processing, particle detection, and filtering stages, and introduces additional analysis logic that operates on the filtered, per-frame position data. This design minimizes duplication, preserves backward compatibility, and allows users to obtain both climbing velocity and fall-related metrics from the same tracked data without modifying upstream components.

# Software Description

## Integration with *FreeClimber*

*FreeClimber-FNG* is implemented in `detector_fng.py`, a standalone fork of the original `detector.py`. The core *FreeClimber* workflow remains unchanged:

1. **Video → cleaned frame stack**
   Background subtraction and cropping.

2. **Frame stack → detected fly positions**
   Particle detection with TrackPy and filtering.

3. **Filtered positions → vial-wise height traces**
   Mean y-position per frame for each vial.

4. **Height traces → FNG events (this extension)**
   Detection of peaks and subsequent descents.

The extension operates entirely on the `df_filtered` table produced in `step_5`, requiring only `frame`, `vial`, and `y` (after FreeClimber's inversion, higher y = higher in vial).

The extension writes a standardized per-event table to `videoname.fng.csv`.

## FNG Detection Logic

For each vial:

1. Compute the mean vertical position for each frame.
2. Apply a centered rolling-mean smoothing kernel.
3. Normalize the trace to a 0–1 range.
4. Detect *climb peaks*.
5. Identify the subsequent fall segment and compute:
   - `frame_peak`
   - `frame_fall_start`
   - `frame_fall_end`
   - `fall_duration_frames`
   - `fall_duration_sec`
   - normalized rise magnitude
   - normalized drop magnitude
   - fall distance in pixels
   - fall distance in centimeters (if pixel scaling available)
   - `recovery_duration_sec` (time from fall end to first sustained climb)

All parameters (smoothing window, minimum rise, minimum drop, minimum separation) are accessible for tuning, enabling adaptation to different rigs and frame rates.

## Implementation

The tool is written in Python and uses the same scientific stack as *FreeClimber*:

- **NumPy**, **pandas** for data manipulation
- **SciPy** (`find_peaks`) for peak detection
- **TrackPy** for object tracking
- **Matplotlib** for optional visualization

The extension does not modify the GUI or underlying detection logic, ensuring backward compatibility. The FNG detection method (_detect_fng_series) can be called on any detector object given a
per-frame vial height trace.

# Quality Control

The FNG detection module was validated on climbing videos (*n* = 120) from laboratory-produced negative geotaxis assays. During development, event detections were cross-checked with the underlying height traces and corresponding video frames to ensure consistency between hand raters and the software.

Thresholds, smoothing windows, and event-gap parameters were tuned until automated detections matched researcher-verified climb–fall behavior. As assay geometry and frame rate vary across laboratories, the software exposes these parameters for user control. Recommended practice is to:

1. Run *FreeClimber* on a small subset of videos.
2. Inspect detected FNG events alongside height traces.
3. Adjust parameters as needed.
4. Apply the final configuration batch-wise across experiments.

To assess agreement between FreeClimber-FNG and manual scoring, a Bland-Altman analysis was conducted on 120 wild-type *Drosophila melanogaster* across three trials each (387 fly-trial observations). Three independent raters scored FNG events from the same videos; hand scores were collapsed across raters and compared to software-detected counts. The mean difference between hand and program scores was approximately zero, with 95% limits of agreement of −1.14 to +0.91 FNG events, indicating that automated detection deviated from manual scoring by less than one event in either direction for the large majority of observations. Discrepancies were concentrated at higher mean event counts, consistent with expected boundary-case disagreement. These results support the use of FreeClimber-FNG as a reproducible alternative to manual scoring in wild-type assays. Validation under additional experimental conditions, including traumatic brain injury and pharmacological manipulation, is ongoing and will be reported in subsequent work.

In addition to human-rater comparison, the pipeline is validated against a synthetic single-fly video with scripted ground-truth fall events (5 events, peak frames known a priori). Automated detection localizes peaks within ±2 frames of ground truth and produces no impossible frame indices. The synthetic video, ground-truth CSV, and regression test are included in the repository under `tests/fixtures/synthetic_validation/`.

![Bland-Altman plot of agreement between hand scoring and FreeClimber-FNG across 387 wild-type fly-trial observations.](freeclimber_fng_bland_altman.png)

# Availability

- **Source Code:** <https://github.com/jordanvasu/FreeClimber-FNG>
- **License:** MIT
- **Operating Systems:** Any system supported by FreeClimber (Windows, macOS, Linux).
- **Dependencies:** Python ≥3.9

This version of the software is archived on Zenodo [@vasu2025freeclimberfng].

# Acknowledgements

This project builds on the original FreeClimber platform by Spierer and colleagues [@spierer2021freeclimber]. I thank Adam Spierer for making the FreeClimber codebase openly accessible. The authors declare no conflicts of interest and received no external funding for this work.

# AI Usage Disclosure

The software was developed and validated by the authors. During development, OpenAI Codex and Anthropic Claude Code were used to audit the codebase for integration issues, logic errors, and security vulnerabilities. All identified issues were reviewed by the authors before fixes were implemented. General-purpose development tools were occasionally used to assist with debugging or refactoring. The algorithm design, parameter selection, validation against experimental data, and integration into the FreeClimber pipeline were performed by the research team.

# References
