<h1>Free Climber</h1>

<h3>Overview</h3>

`FreeClimber` was originally designed to quantify the clmimbing velocity of fruit flies in a negative geotaxis assay. This program can also be co-opted for subtracting the static background of a video, as well as particle detection of other videos where points of uniform size on a lighter background move from the bottom to top of a video. 

This gif demonstrates the principle underlying how velocity is quantified:

<img src="https://github.com/adamspierer/FreeClimber/blob/master/z/0-Tutorial_climbing.gif" width="500" height="167" align="center">

On the left, blue '+' represent candidate spots, while colored circles represent 'true' spots that passed a specified threshold. Circles are blue during the most linear portion of the position vs. time (velocity) curve, and red otherwise. On the right, the mean vertical position of all points by frame show a gray line connecting the mean vertical position for all true spots identified by frame. The line of best fit corresponds with the regression line through the most linear segment, and its slope corresponds with the climbing velocity of that cohort.

This platform implements a Graphical User Interface (GUI) so users can optimize detection parameters and specify experimental details, and a command line interface for batch processing many videos with the same detection parameters.

<h2>Installation and requirements</h2>

<h3>Requirements</h3>

    - argparse      [1.1   ]
    - ffmpeg-python [0.2.0 ]
    - matplotlib    [3.1.3 ]
    - numpy         [1.18.1]
    - pandas        [1.0.0 ]
    - pip           [20.0.2]
    - scipy         [1.4.1 ]
    - trackpy       [0.4.2 ]
    - wxPython      [4.0.4 ]

NOTE: We recommend using with a Python3.6 virtual environment, though it was built with a Python3.7 environment.

<h3>Installing</h3>

We recommend running this package in an Anaconda-based virtual environment. Anaconda can be downloaded [here](https://docs.anaconda.com/anaconda/install/).


**Make sure `conda` is installed** (should return something like `conda 4.7.11`):

	conda -V 

**Update conda if needed** (press `y` when prompted):

	conda update conda

**Create a Python 3 virtual environment** (replace `python36` with your name of choice):
	
	conda create -n python36 python=3.6 anaconda

*NOTE: See note above about Python 3.6 vs. 3.7*

**Activate your virtual environment**:

	conda activate python36
	
**OR** (if that doesn't work):

	source activate python36

For more details about creating a conda virtual environment, see [here](https://uoa-eresearch.github.io/eresearch-cookbook/recipe/2014/11/20/conda/). Once the environment is set up and activated, we can install the dependencies listed in the `Requirements` section above.

**Installation using PyPi** <!-- (Recommended)**: -->

	pip install FreeClimber
	pip install ffmpeg-python

<!--
**Using conda (still in testing)**:

    conda install -c adamspierer freeclimber
-->

**Download the script files** (can be done with `git clone` if user is familiar with `git` or by directly downloading the `.py` files into a single folder)

**Cloning the git repository**:

	cd <folder of interest>
	git clone https://github.com/adamspierer/FreeClimber.git
	
NOTE: As of now, the platform itself is <u>not</u> callable as a module and these steps merely download the dependencies. The script files must be directly referenced when running the program. See our [tutorial](https://github.com/adamspierer/FreeClimber/blob/master/TUTORIAL.md) for usage instructions.

<h2>Running the interactive, Graphical User Interface</h2>

The graphical user interface (GUI) is useful for optimizing the parameters one might use in analyzing a video. Several fields are available to modify and buttons enable quick testing.

Before getting started, it is important that files are named with the experimental variables in the file name and separated by underscores. It is also important that variables attached to a specific vial (`Vial_ID variables`) are in the beginning. Here, genotype and sex are unique to each vial, but the day and replicate change over the experiment.
	
	Example: <genotype>_<sex>_<day>_<replicate>.<suffix> ==> w1118_m_2.h264


Additionally, the adage: "garbage in, garbage out" applies to these videos. Higher quality videos with clean backgrounds and strong contrast between spots and the background will yield the best results. FreeClimber can work with other videos, but they may result in noiser results.

<h3>Step 1a - Specify a video </h3>

To begin, we will start from the FreeClimber folder and run the GUI to optimize the detection parameters: 

	pythonw ./scripts/FreeClimber_gui.py --video_file ./example/w1118_m_2_1.h264

If you get an error like the following, then be sure to call `pythonw` instead of `python`:
    
    This program needs access to the screen. Please run with a Framework build of python, and only when you are logged in on the main display of your Mac.

If you get an error like the following then be sure to run the script in the (correct) virtual environment:
	
	ModuleNotFoundError: No module named 'ffmpeg'

<img src="https://github.com/adamspierer/FreeClimber/raw/master/z/1-opening.png" width="500" height="400">

Alternatively, you can omit the `--video_file` flag and path for a dialog box to appear for browsing files. The imported file path will be in a box between the control panel and visualization window.

You can also select the `Browse...` button in the upper left to chose a different video at any time.


<h3>Step 1b - Define options </h3>

When the video loads, the program will collect what metadata it can and populate various fields across the control panel.

In this step, we are specifying the "Pixels / cm" (`pixels_to_cm` variable) and the "Frames / sec" (`frame_rate` variable).

**Pixels / cm** - Calculate by drawing a horizontal or vertical line between two points of known distance. In this example, we can a horizontal line to measure the width of a [narrow glass vial](https://geneseesci.com/shop-online/product-details/32-201/drosophila-vials-narrow-glass2?search=glass%20vials) (2.45 cm). Record the value in the `Width` or `Height` field for each of the vials, calculate the average, and divide that by the known distance in cm. 

<img src="https://github.com/adamspierer/FreeClimber/raw/master/z/2a-measuring.png" width="500" height="400">

Here, we measured three vials at 110, 114, 111 pixels (average = 111.7 pixels/vial). When we divide that by the known width of a vial (2.45 cm), we are left with a conversion factor of 45.6 pixels/cm.

<img src="https://github.com/adamspierer/FreeClimber/raw/master/z/2b-pixel_to_cm_calculation.png" width="223" height="200">

**Frames / sec** - The default frame rate is 25, but may not be the actual frame rate. If you are unsure of the recording device's frame rate, record a video over a fixed time, load that video, and divide the second field in the `Background Frames` (number of frames) by the recording time to get the frame rate.

If you decide to do these steps, make sure to click the `Convert to cm / sec` box so final slopes are converted to `cm / second`. Otherwise, the resulting slopes will be in `pixels / frame`.

<h3>Step 2 - Select ROI </h3>

The following fields are populated with the video metadata and result in the entire image being analyzed. You can click and drag on the image from the upper-right to left-left to specify a smaller Region of Interest (ROI). This will drastically improve performance.

**x-pos.** - Left-most point to crop the Region of Interest (ROI). Default is 0.

**y-pos.** - Top-most point to crop the ROI. Default is 0.

**Width** - The number of pixels in the `x` direction. Starting value is the width of the video, but this changes when selecting a region of interest. Default is image width dimension.

**Height** - The number of pixels in the `y` direction. Starting value is the height of the video, but this changes when selecting a region of interest. Default is image height dimension.

<img src="https://github.com/adamspierer/FreeClimber/raw/master/z/3-boxing.png" width="500" height="400">

A checkbox for **Fixed ROI size?** allows the user to define specific dimensions. Once the dimensions are specified and the box is checked, users can click on the image to move the box around.

NOTE: Images are indexed from upper left to lower right, which is different from most plots which are oriented lower left to upper right. Think keypads on a phone vs. calculator.

<h3>Step 3 - Specify spot parameters</h3>

These parameters configure the detector and filter for candidate vs. true spots based on a `minmass`, signal `threshold`, and spot eccentricity/`ecc`/circularity. If you are processing many videos, it is best practice to optimizing these parameters with a few representative videos and check out the `.spot_check.png` file to determine if the final three variables are suitable.

**Diameter** - Approximate diameter for estimating a spot. Only odd numbers, rounding up when in doubt.

**MaxDiameter** - The maximum distance across to consider a spot.

**MinMass** - The minimum integrated brightness. Spots will be filtered if their mass value does not pass this threshold.

**Threshold** - The signal threshold value for filtering spots. Specify a threshold, or leave as "auto" for the program to calculate one. The "auto" threshold looks for the local minimum between two local maximum, or takes the halfway mark between the global maximum.

**Ecc/circularity** - Each spot has a measure of its circularity (0 is more circular). In our experience, the upper limit is more important, and True spots don't often exceeded 0.6.

For more information, see the [trackpy documentation](http://soft-matter.github.io/trackpy/v0.4.2/generated/trackpy.locate.html#trackpy.locate).

<h3>Step 4 - Additional parameters</h3>

These parameters are important for setting the `background frames` range for calculating the background image, `checking frames` during the optimization process, specifying the number of `vials` in a video, and the `window size` of the local linear regression. 

**Background Frames** - The first through last frames from which to calculate a background image. The default is to use all frames, but a subset can be used to account for unexpected background movement or irregularities in the video. Here, the median pixel intensity for each x,y-coordinate is calculated. This median image (middle-left image) is subtracted from each individual frame (middle), resulting in videos with much greater signal:noise ratio (middle-right image; right image). This background-subtracted video can then be processed by the detector (left-most image).

<img src="https://github.com/adamspierer/FreeClimber/raw/master/z/4-background_subtraction.png" width="600" height="200">


**Crop Frames** - The first and last frames to include if the video is too long. `Background frames` and `Check frames` should start from the first frame loaded, not the first frame after cropping.

**Check Frame** - Enter a specific frame to display the candidate (semi-transparent blue +) and True points (colored circles) when running the `Test Parameters` button. The default is the first frame. 

NOTE: This field is only important for the GUI.

**Number of vials** - The number of vials present in a video. Spots are assigned to lanes or vials based on their x-position. These divisions are evenly-spaced bins placed between the maximum and minimum x-values for all frames. Some issues can arise if there are true spots (crossing signal and mass threshold) that lie beyond the range that they should (ex. true spot as an artefact of the rig) or if in a batch of videos one vial has no flies that move (ex. longitudinal aging study). To circumvent these issues, you can make `Number of vials = 1`.

**Window size** - Number of frames for a sliding window when calculating the local linear regression. The best practice is to input ~90% of the number of frames corresponding with the most linear portion of the asymptotic or sigmoidal mean y-position vs. time plot for one of the fastest climbers in your panel. In our work, we found flies can climb the vial in about seconds at 29 frames per second `29 * 2 * 0.9 = 52`.

**Trim outliers?** - Check if you want to do a smart trim of the data. This is useful if there are outliers on the periphery. Values for the Top and Bottom (**TB**) and Left and Right (**LR**) correspond with sensitivity scores for determining if a point or groups of points on the extreme edges of the image are outliers or not. Leave the box unchecked to ignore.

<h3>Step 5 - Naming parameters</h3>

These next variables are important for properly saving experimental details and file paths.

**Naming pattern** - As noted above, all input files MUST be named with experimental details separated by an underscore. The program will scrap these details from the file name and fill them into resulting data file columns named by the `Naming pattern`. 

If an experiment included testing several genotypes of both sexes across days with several technical replicates for each, then the following would be appropriate:

	Naming pattern = genotype_sex_day_rep

Experimental details that do not change over the course of the experiment (e.g. genotype and sex) should be in the beginning with the number of these fields entered in **Vial ID variables**. Experimental details that change (e.g. time and replicate) should be at the end. Within each of these distinctions, the order does not matter.

**Vial_ID variables** - The number of variables that are consistent during the experiment (e.g. 2 for genotype and sex). These become the variables associated with an individual vial over a longitudinal study. See entry above.

**Project path** - The parent folder to search through for all video files of a given type (same video file type as the video file loaded). The generated configuration (`.cfg`) file and final `results.csv` file containing all slopes from all videos will eventually be saved to this folder.

<h3>Buttons</h3>

When all the appropriate fields are set, we can process the video.

**Test parameters** - This button begins a full analysis of the loaded video and plots out several diagnostic plots. On the top row (top row, left to right) the median background image, `Check frame` frame number with candidate (blue +; none here) and true (colored circles) spots, and the mean vertical-position vs. time plots (darker segments indicate most linear section). On the bottom row (left to right), there are plots for the distribution of spots' mass (line =  `minmass` value), distribution of spot signals (line = signal `threshold` value), and number of spots counted per frame (for entire video). 

To make adjustments, modify the appropriate field(s) and `Test parameters`, or select a new video from `Browse...`. If the program freezes due to poor spot quality, press the `Reload video` button.

When all variables are appropriately filled, press `Save parameters` to generate the final `.cfg` file.

<img src="https://github.com/adamspierer/FreeClimber/raw/master/z/5-display.png" width="500" height="400">

**Reload video** - Reload a video to refine detection parameters once the `Test parameters` button has been pressed, or load a new video with the **Browse** button. This button is also helpful if the program hangs up.

**Save parameters** - Generate a `.cfg` file from the displayed parameters. This file can be edited manually to use the same analysis parameters on a different set of files. Just make sure to press this when you are done optimizing, otherwise no file will save.

<h3>Optimization plots</h3>

To facilitate researchers' ability to optimize `FreeClimber`, the GUI will also generate three figure files with the following suffixes:

1. `ROI.png` - First frame of the color video with the region of interest outlined in a red rectangle. The white lines down the middle represent the predicted placement of the vial dividers, based on the distribution of flies within the frame for a more informed estimate of vial placement. The cyan lines represent the ROI after adjusting for trimming outliers, if it is selected.
    
<img src="https://github.com/adamspierer/FreeClimber/blob/master/example/w1118_m_2_1.ROI.png" width="500" height="400">    
    
2. `processed.png` - Designed to help with a frame range for the `Background frames`. Three plots demonstrate what the video looks like cropped and grayscaled, as a null background, and with the background subtracted. The frame number corresponds with the "Check frame"/`check_frame` field/variable.
    
<img src="https://github.com/adamspierer/FreeClimber/blob/master/example/w1118_m_2_1.processed.png" width="375" height="500">

3. `spot_check.png` - Three sets of two subplots, corresponding with the distribution of certain spot metrics and where the spots lie on the image. This is especially useful for determining eccentricity, mass, and signal thresholds.
    
<img src="https://github.com/adamspierer/FreeClimber/blob/master/example/w1118_m_2_1.spot_check.png" width="522" height="300">  

<h2>Running the command line version</h2>

To use, run: 

`python ./scripts/FreeClimber_main.py --config_file ./example/example.cfg`

This command differs from the GUI in a few ways. First, it uses `python` instead of `pythonw`. Second, it has a `--config_file` flag instead and looks for a file with the `.cfg` suffix. 

Additionally, using the command line interface is useful for batch processing subsets of files. The default `--process_all` will process all files in the `path_project` with a common `file_suffix`. If new files are added, you can call `--process_undone` to process files that don't have the corresponding `.slopes.csv` (individual video) results file. 

If you only want to process certain files, use the `--process_custom` flag followed by the path to a processed (`.prc`) file:

`python ./scripts/FreeClimber_main.py --config_file ./example/example.cfg --process_custom ./example_other/custom.prc`

The `.prc` file is a text file with the file paths to videos of interest. One can be generated by the user, so long as it has a `.prc` suffix. Alternatively, the `gather_files.py` script will create a file with all appropriate file paths in the parent folder that can be cleaned up in a text editor:

`python ./scripts/gather_files.py --parent_folder ./example_other --suffix h264 --save_file`

We also provide flags for `--optimization_plots` (generates files with the `spot_check.png`, `ROI.png`, and `processed.png` suffixes for optimizing the detection parameters, region of interest, and background subtraction parameters, respectively). Though this will do so for every video when run through the command line.

For each of the scripts provided, help documentation is provided if you type:

    python <path_to_file.py> -h

<h4>Additional test videos</h4>

There are five additional test videos in the `example_other` folder. These videos have two from the same rig set up as the main test/example. One of these has all six vials filled (climbing\_1), while the other is similar but demonstrates how the program works well when there is a wagon-wheel/banding effect (climbing\_2). 

The remaining three videos represent how an improvised rig might perform, in order of decreasing quality. In the first of these, (climbing\_3) has a clean background and uses a custom rig (piece of wood with even spaced vial-diameter holes drilled in the top, held together with a hinge and rubberbands). This is the same video from the gif at the beginning of the README and TUTORIAL files. The next video (climbing\_4) uses a test tube rack as the vial holder, has a dramatic shift of the rig in the first second, and includes partially erased sharpie marks on the vials. These marks might otherwise be false positives without the background subtraction step and movement in the first frames requires that the video is cropped in the program before processing. The third of these videos (climbing\_5) is on the lower end for quality in what FreeClimber can handle. It has four vials (three with flies) strapped to a piece of plywood and partially resting on a cardboard box. The vials have some abrasion so there is not always a clear shot of the flies. While this video is far from an ideal input, it is nonetheless impressive that data can be extracted from the video.

<h4>Tips</h4>

- We recommend trimming the length of the video down to the length you are interested in analyzing and having a single replicate of the assay per video. To slice videos efficiently using a command line-based tool, we recommend using `FFmpeg`, a well documented media manipulation tool: https://trac.ffmpeg.org/wiki/Seeking. Alternatively, frames to trim can be specified in the `Crop frames` fields of the GUI.

- For best spot-detection, we also recommend videos have a clean and consistent background. One of the advantages to this platform is the method of background subtraction. As long as the background is static and there is sufficient contrast between it and the flies, the detector should work reasonably well. We found backlighting the flies with an LED-tracing board created strong contrast and made for better data quality. 

- If using an LED-tracking board, one issue that can arise is a 'wagon-wheel effect' where a "rolling shutter" paired with a repetitive/cyclical phenomenon creates an optical illusion. In this case, the LEDs will create a horizontal banding effect that will have a slight negative affect the ability of the detector to subtract the background and negatively affect the data quality. Adjust the dimmer settings until you find a brightness that works.

- When optimizing spot detection parameters, first focus on the `Diameter`, `MaxDiameter`, and `MinMass`. These fields determine whether a group of pixels is a candidate spot so find the highest values for these parameters that cover all obvious spots, since lower values can increase the time required for the program to process a video. The `Ecc` field should be next, followed by the `Threshold` field. When `Threshold` is set to "auto" it will either find the local minima between two peaks, or divide the global peak score in half. The `<file_name.spot_check.png>` wil be helpful for optimizing the detector parameters. 

- Image coordinates are numbered from upper left to lower right, so <file_name.raw.csv> and `<file_name.filtered.csv>` will have y-coordinates inverted for plotting purposes. This is due to how images are indexed and graphs plotted, and is important for those looking to use the outputs for further visualization.

- The `<file_name.raw.csv>` output is also really useful for passing to `trackpy` for linking spots together into tracks. Read more on TrackPy particle __tracking__ [API](http://soft-matter.github.io/trackpy/v0.4.2/api.html). Briefly, the file can be loaded in and saved as `f`, which will seamlessly integrate into the `Step 3: Link features into particle trajectories` in the [TrackPy Walkthrough](http://soft-matter.github.io/trackpy/v0.4.2/tutorial/walkthrough.html).

<h4>Table of variables</h4>

|Variable name |	Data type |	Explanation|
| --- | --- | --- |
|x |	Integer* |	Leftmost pixel of ROI for analysis|
|y |	Integer*	| Topmost pixel of ROI for analysis|
|w | Integer* |	Width of ROI for analysis|
|h | Integer* |	Height of ROI for analysis|
|check\_frame |	Integer	| First frame to display in Test Parameters|
|blank\_0 |	Integer	| First frame of range to subtract background|
|blank\_n |	Integer	| Last frame of range to subtract background|
|crop\_0 |	Integer	| First frame of range to crop|
|crop\_n |	Integer	| Last frame of range to crop|
|threshold |	Integer	| Threshold for filtering against points|
|diameter	|Integer	| Estimated diameter of spot in pixels|
|minmass	|Integer	| The minimum integrated brightness|
|maxsize | Integer| The maximum distance across to consider a spot|
|ecc\_low |	Float	| Lower bounds for spot circularity (0 - 1) |
|ecc\_high | Float	| Upper bounds for spot circularity (0 - 1) |
|vials	| Integer |	Number of vials in video|
|window |	Integer	| Number of frames for sliding window|
|pixel\_to\_cm|	Integer*	| Conversion factor for pixels to centimeters|
|frame\_rate	|Integer* |	Video frame rate|
|vial\_ID\_vars	|Integer |	Number of variables in naming convention that are consistent across a time-dependent experiment (ex. genotype, sex)|
|outlier\_TB |	Integer* |	Top and bottom sensitivity factor for trimming outliers|
|outlier\_LR |	Integer* |	Left and right sensitivity factor for trimming outliers|
|naming\_convention |	String |	Experimental conditions in file name|
|path\_project| String |	Path to parent folder containing experimental files, configuration\_file.cfg and end results.csv file eventually saved here|
|file\_suffix |	String |	Suffix of videos being processed|
|convert\_to\_cm\_sec |	Boolean | True if converting output slope to centimeters per second|
|trim\_outliers |	Boolean | True if trimming outliers|

\* - Can be either an integer or float