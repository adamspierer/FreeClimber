<h1>Free Climber</h1>

<h3>Overview</h3>

`FreeClimber` was originally designed to quantify the clmimbing velocity of fruit flies in a negative geotaxis assay. This program can also be co-opted for subtracting the static background of a video, as well as particle detection of other videos where points of uniform size on a lighter background move from the bottom to top of a video. 

This gif demonstrates the principle underlying how velocity is quantified:

<img src="https://github.com/adamspierer/FreeClimber/blob/master/z/0-Tutorial_climbing.gif" width="500" height="167" align="center">

On the left, blue '+' represent candidate spots, while colored circles represent 'true' spots that passed a specified threshold. Circles are blue during the most linear portion of the position vs. time (velocity) curve, and red otherwise. On the right, the mean vertical position of all points by frame show a gray line connecting the mean vertical position for all true spots identified by frame. The line of best fit corresponds with the regression line through the most linear segment, and its slope corresponds with the climbing velocity of that cohort.

This platform implements a Graphical User Interface (GUI) so users can optimize detection parameters and specify experimental details, and a command line interface for batch processing many videos with the same detection parameters.

<h2>Installation and requirements</h2>

We recommend running this package in an Anaconda-based virtual environment. Anaconda can be downloaded [here](https://docs.anaconda.com/anaconda/install/).

**Make sure `conda` is installed** (should return something like `conda 4.7.11`):

	conda -V 

**Update conda if needed** (press `y` when prompted):

	conda update conda

**Create a Python 3 virtual environment** (replace `python36` with your name of choice):
	
	conda create -n python36 python=3.6 anaconda

*NOTE: See note above about Python 3.6 vs. 3.7..*

**Activate your virtual environment**:

	conda activate python36
	
**OR** (if that doesn't work):

	source activate python36

For more details about creating a conda virtual environment, see [here](https://uoa-eresearch.github.io/eresearch-cookbook/recipe/2014/11/20/conda/). Once the environment is set up and activated, we can install the dependencies listed in the `Requirements` section above.

Once the environment is set up and activated, we can install the following dependencies.

    - argparse   [1.1   ]
    - ffmpeg     [1.4   ]
    - matplotlib [3.1.3 ]
    - numpy      [1.18.1]
    - pandas     [1.0.0 ]
    - pip        [20.0.2]
    - scipy      [1.4.1 ]
    - trackpy    [0.4.2 ]
    - wxPython   [4.0.4 ]

NOTE: This platform was developed using Python 3.7.6, and tested in a Python 3.6 environment. We recommend running in a Python 3.6 environment for best results, earlier Python versions are untested.

**Installation using PyPi (Recommended)**: 

	pip install FreeClimber
    pip install ffmpeg-python

NOTE: We recommend downloading these two modules separately.

**Installation using `conda` (still in beta-testing)**:

	conda install -c adamspierer freeclimber


**Download the script files** (can be done with `git clone` if user is familiar with `git` or by directly downloading the `.py` files from the GitHub repository and placed into a single folder)

**Cloning the git repository**:

	cd <folder of interest>
	git clone https://github.com/adamspierer/FreeClimber.git
	
NOTE: As of now, the platform itself is <u>not</u> a module and these steps merely download the dependencies. The script files must be directly referenced when running the program. The files in the Github repository are the most up-to-date, while those connected to the PyPi or conda distributions are not.

<h2>Running the interactive, Graphical User Interface</h2>

The graphical user interface (GUI) is useful for optimizing the parameters one might use in analyzing a video. Several fields are available to modify and buttons enable quick testing.

Before getting started, it is important that the videos have spots with sufficient contrast against a light background, and that files are named with the experimental variables in the file name and separated by underscores. 
	
	Example: <genotype>_<sex>_<replicate>.<suffix> ==> w1118_m_2.h264

<h3>Step 1 - Specify a video</h3>

To begin, we will run the GUI to optimize the detection parameters: 

	pythonw ./scripts/FreeClimber_gui.py --video_file ./example/w1118_m_2_1.h264

If you get an error like the following, then be sure to call `pythonw` instead of `python`:
    
    This program needs access to the screen. Please run with a Framework build of python, and only when you are logged in on the main display of your Mac.


<img src="https://github.com/adamspierer/FreeClimber/raw/master/z/1-opening.png" width="500" height="400">

Alternatively, you can omit the `--video_file` flag and path for a dialog box to appear for browsing files. The imported file path will be in a box between the control panel and visualization window.

When the video loads, the program will collect what metadata it can and populate the fields in the control panel:

**x-pos.** - Left-most point to crop the Region of Interest (ROI).

**y-pos.** - Top-most point to crop the ROI.

**Width** - The number of pixels in the `x` direction. Starting value is the width of the video, but this changes when selecting a region of interest

**Height** - The number of pixels in the `y` direction. Starting value is the height of the video, but this changes when selecting a region of interest

NOTE: Images are indexed from upper left to lower right, which is different from most plots which are oriented lower left to upper right. The loaded video will show the first frame, so no need to be alarmed when the y-axis reflects this aspect of image indexing.

**Frame rate** - The number of `frames / second` has a value of 25. This value is useful if you want to convert the end slope to `cm / second` rather than `pixels / frame`.

<h3>Step 2 - Defining video parameters</h3>

These variables are important for adding a scaling factor if you are interested in having `cm / seconds` vs `pixels / frame` (advantage being that results can be compared across studies) and decreasing the time required to analyze a single video.

To start with, we can select-and-drag the cursor across the video to draw a rectangle. This rectangle can be redrawn several times, making it useful for:
	
1. <i><u>Measuring known distances to calibrate the final slopes into pixels per centimeter</i></u> 

	To calculate the **Pixels per centimeter** parameter, draw a horizontal or vertical line between two points of known distance. In this example, we can a horizontal line to measure the width of a [narrow glass vial](https://geneseesci.com/shop-online/product-details/32-201/drosophila-vials-narrow-glass2?search=glass%20vials) (2.45 cm). Record the value in the `Width` or `Height` field for each of the vials, calculate the average, and divide that by the known distance in cm. 

	<img src="https://github.com/adamspierer/FreeClimber/raw/master/z/2a-measuring.png" width="500" height="400">

	Here, we measured three vials at 110, 114, 111 pixels (average = 111.7 pixels/vial). When we divide that by the known width of a vial (2.45 cm), we are left with a conversion factor of 45.6 pixels/cm. For our example, we will use 45.

	<img src="https://github.com/adamspierer/FreeClimber/raw/master/z/2b-pixel_to_cm_calculation.png" width="223" height="200">

	To automatically convert the final slopes to `cm / second`, click the `Convert to cm / sec` box. Alternatively, the resulting slopes will be in `pixels / frame`.

2. <u><i>Drawing a Region of Interest (ROI) to improve detection time.</i></u>

	Next, we can specify a ROI to improve the performance of the detector by only analyzing the parts of the video we care about. Draw a rectangle from upper-left to lower-right. The text box at the bottom will remind you, but other directions will lead to negative values in the positional fields will result in an error.

	<img src="https://github.com/adamspierer/FreeClimber/raw/master/z/3-boxing.png" width="500" height="400">

	A checkbox for **Fixed ROI size?** allows the user to define specific dimensions. Once the dimensions are specified and the box is checked, users can click on the image to move the box around.

<h3>Step 3 - Specify spot parameters</h3>

These parameters are for configuring the detector and filtering for candidate vs. true spots based on a `minmass` and signal `threshold`. If you are processing many videos, it is best practice to optimizing these parameters with a few representative videos.

**Diameter** - Approximate diameter for estimating a spot. Only odd numbers, rounding up when in doubt.

**MaxDiameter** - The maximum distance across to consider a spot.

**MinMass** - The minimum integrated brightness. Spots will be filtered if their mass value does not pass this threshold.

**Threshold** - The signal threshold value for filtering spots. Specify a threshold, or leave as "auto" for the program to calculate one. The "auto" threshold looks for the local minimum between two local maximum, or takes the halfway mark between the global maximum.

For more information, see the [TrackPy documentation](http://soft-matter.github.io/trackpy/v0.4.2/generated/trackpy.locate.html#trackpy.locate).

<h3>Step 4 - Additional parameters</h3>

These parameters are important for setting the `background frames` range for calculating the background image, `checking frames` during the optimization process, specifying the number of `vials` in a video, and the `window size` of the local linear regression. 

**Background Frames** - The first through last frames from which to calculate a background image. The default is to use all frames, but a subset can be used to account for unexpected background movement or irregularities in the video. Here, the median pixel intensity for each x,y-coordinate is calculated. This median image (middle-left image) is subtracted from each individual frame (middle), resulting in videos with much greater signal:noise ratio (middle-right image; right image). This background-subtracted video can then be processed by the detector (left-most image).

<img src="https://github.com/adamspierer/FreeClimber/raw/master/z/4-background_subtraction.png" width="600" height="200">

**Check Frames** - Enter a specific frames to display the candidate (yellow x) and True points (colored circles) when running the `Test Parameters` button. The defaults is the first frame. 

NOTE: This field is only important for the GUI.

**Number of vials** - The number of vials present in a video. Spots are assigned to lanes or vials based on their x-position. These divisions are evenly-spaced bins placed between the maximum and minimum x-values for all frames. Some issues can arise if there are true spots (crossing signal and mass threshold) that lie beyond the range that they should (ex. true spot as an artefact of the rig) or if in a batch of videos one vial has no flies that move (ex. longitudinal aging study). To circumvent these issues, you can make `Number of vials = 1`.

**Window size** - Number of frames for a sliding window when calculating the local linear regression. The best practice is to input ~90% of the number of frames corresponding with the most linear portion of the asymptotic or sigmoidal mean y-position vs. time plot for one of the fastest climbers in your panel. In our work, we found flies can climb the vial in about seconds at 29 frames per second `29 * 2 * 0.9 = 52`, so we settled on 50 frames.

<h3>Step 5 - Naming parameters</h3>

These next variables are important for properly saving experimental details and file paths.

**Naming pattern** - As noted above, all input files should be named with experimental details separated by an underscore. The program will scrap these details and fill them into columns named by the `Naming pattern` in the appropriate results files. For instance, if an experiment was testing several genotypes of both sexes in a longitudinal study across days with several technical replicates for each, then:

	Naming pattern = genotype_sex_day_rep
	
The program will automatically load the file name, just be sure to replace these fields with the appropriate variable name. There is no need for a file suffix since this is pulled from the input file and the order is somewhat important. Experimental details that do not change over the course of the experiment (e.g. genotype and sex) should be in the beginning with the number of these fields entered in **Vial ID variables**. Experimental details that change (e.g. time and replicate) should be at the end. Within each of these distinctions, the order does not matter.

**Vial_ID variables** - The number of variables that are consistent during the experiment (e.g. 2 for genotype and sex). See entry above.

**Project path** - The parent folder to search through for all video files of a given type (same video file type as the video file loaded). The generated configuration (`.cfg`) file and final `results.csv` file containing all slopes from all videos will eventually be saved to this folder.

<h3>Buttons</h3>

When all the appropriate fields are set, we can process the video.

**Test parameters** - This button begins a full analysis of the loaded video and plots out several diagnostic plots. On the top row (top row, left to right) the median background image, `Check frame` frame number with candidate (yellow x; none here) and true (colored circles) spots, and the mean vertical-position vs. time plots (darker segments indicate most linear section). On the bottom row (left to right), there are plots for the distribution of spots' mass (line =  `minmass` value), distribution of spot signals (line = signal `threshold` value), and number of spots counted per frame (for entire video). 

To make adjustments, modify the appropriate field(s) and `Test parameters`, or select a new video from `Browse...`. If the program freezes due to poor spot quality, press the `Reload video` button.

When all variables are appropriately filled, press `Save parameters` to generate the final `.cfg` file.

<img src="https://github.com/adamspierer/FreeClimber/raw/master/z/5-display.png" width="500" height="400">

**Reload video** - Reload a video to refine detection parameters once the `Test parameters` button has been pressed, or load a new video with the **Browse** button.

**Save parameters** - Generate a `.cfg` file from the displayed parameters. This file can be edited manually to use the same analysis parameters on a different set of files.

<h3>Optimization plots</h3>

To facilitate researchers' ability to optimize `FreeClimber`, the GUI will also generate three figure files:

    1. `<file_name>.ROI.png` - First frame of the color video with the region of interest outlined in a red rectangle. The white lines down the middle represent the predicted placement of the vial dividers, based on the distribution of flies within the frame for a more informed estimate of vial placement.
    
<img src="https://github.com/adamspierer/FreeClimber/blob/master/example/w1118_m_2_1.ROI.png" width="500" height="400">    
    
    2. `<file_name>.processed.png` - Designed to help with a frame range for the `Background frames`. Three plots demonstrate what the video looks like cropped and grayscaled, as a null background, and with the background subtracted.
    
<img src="https://github.com/adamspierer/FreeClimber/blob/master/example/w1118_m_2_1.processed.png" width="500" height="400">

    3. `<file_name>.spot_check.png` - Three sets of two subplots, corresponding with the distribution of certain spot metrics and where the spots lie on the image. This is especially useful for determining signal and mass thresholds.
    
<img src="https://github.com/adamspierer/FreeClimber/blob/master/example/w1118_m_2_1.spot_check.png" width="500" height="400">  

<h2>Running the command line version</h2>

To use, run: 

```python ./scripts/FreeClimber_main.py --config_file ./example/example.cfg```

This command differs from the GUI in a few ways. First, it uses `python` instead of `pythonw`. Second, it has a `--config_file` flag instead and looks for a file with the `.cfg` suffix. 

Additionally, using the command line interface is useful for batch processing subsets of files. The default `--process_all` will process all files in the `path_project` with a common `file_suffix`. If the program exits (not because of an erro) part way through a large job, you can call `--process_undone` to process files that don't have the corresponding `.slopes.csv` (individual video) results file. 

If the program exits because of an unknown error or you want to process a select group of videos, use the `--process_custom` flag followed by the path to a processed (`.prc`) file:

    python ./scripts/FreeClimber_main.py --config_file ./example/example.cfg --process_custom ./example_other/custom.prc

This is only an example command and the `custom.prc` file is only an example file, but both illustrate the basic idea. The file is simply a text file with the file paths to videos of interest. One can be generated by the user, so long as it has a `.prc` suffix. Alternatively, this one was generated with the `gather_files.py` script using the following command, and then manually removing the files we don't care for in a text editor:

    python ./scripts/gather_files.py --parent_folder ./example_other --suffix h264 --save_file

We also provide flags for `--optimization_plots` (generates files with the `spot_check.png`, `ROI.png`, and `processed.png` suffixes for optimizing the detection parameters, region of interest, and background subtraction parameters, respectively). Though this will do so for every plot.

For each of the scripts provided, help documentation is provided if you type:

    python <path_to_file.py> -h

<h4>Tips</h4>

- We recommend trimming the length of the video down to the length you are interested in analyzing and having a single replicate of the assay per video. To slice videos efficiently using a command line-based tool, we recommend using `FFmpeg`, a well documented media manipulation tool: https://trac.ffmpeg.org/wiki/Seeking.

- For best spot-detection, we also recommend videos have a clean and consistent background. One of the advantages to this platform is the method of background subtraction. As long as the background is static and there is sufficient contrast between it and the flies, the detector should work reasonably well (see the `noisy_background_few_flies.h264` example, which will require playing with the controls (ex. x = 500, y = 360, Width = 235, Height = 174, Diameter = 5, MaxDiameter = 9, MinMass = 30, Background frames = 70 and 150, all other variables are okay to be left). This video is on the lower end of quality one can hope to get results from. We found backlighting the flies with an LED-tracing board created strong contrast and made for better data quality. 

- If using an LED-tracking board, one issue that can arise is a 'wagon-wheel effect' where a "rolling shutter" paired with a repetitive/cyclical phenomenon creates an optical illusion. In this case, the LEDs will create a horizontal banding effect that will affect the ability of the detector to subtract the background and negatively affect the data quality. Adjust the dimmer settings until you find a brightness that works. See `./examples_other/clean_background_many_flies_with_lines.h264` for a example.

- When optimizing spot detection parameters, first focus on the `Diameter`, `MaxDiameter`, and `MinMass`. These fields determine whether a group of pixels is a candidate spot so find the highest values for these parameters that cover all obvious spots, since lower values can increase the time required for the program to process a video. The `Threshold` field should be optimized last or set to "auto" to determine which spots pass a signal threshold. The <file_name.spot_check.png> wil be helpful for optimizing the detector parameters. 

- Image coordinates are numbered from upper left to lower right, so <file_name.raw.csv> and <file_name.filter.csv> will have y-coordinates inverted for plotting purposes. This is due to how images are indexed and graphs plotted, and is important for those looking to use the outputs for further visualization.

- The `<file_name.raw.csv>` output is also really useful for passing to `TrackPy` for linking spots together into tracks. Read more on TrackPy particle __tracking__ [API](http://soft-matter.github.io/trackpy/v0.4.2/api.html). Briefly, the file can be loaded in and saved as `f`, which will seamlessly integrate into the `Step 3: Link features into particle trajectories` in the [TrackPy Walkthrough](http://soft-matter.github.io/trackpy/v0.4.2/tutorial/walkthrough.html).

<h4>Table of variables</h4>
|Variable name |	Data type |	Explanation|
| --- | --- | --- |
|x |	Integer |	Leftmost pixel of ROI for analysis|
|y |	Integer	| Topmost pixel of ROI for analysis|
|w | Integer |	Width of ROI for analysis|
|h | Integer |	Height of ROI for analysis|
|frame_0 |	Integer	| First frame to display in Test Parameters|
|blank_0 |	Integer	| First frame of range to subtract background|
|blank_n |	Integer	| Last frame of range to subtract background|
|threshold |	Integer	| Threshold for filtering against points|
|diameter	|Integer	| Estimated diameter of spot in pixels|
|minmass	|Integer	| The minimum integrated brightness|
|maxsize | Integer| The maximum distance across to consider a spot|
|vials	| Integer |	Number of vials in video|
|window |	Integer	| Number of frames for sliding window|
|pixel_to_cm|	Integer	| Conversion factor for pixels to centimeters|
|frame_rate	|Integer |	Video frame rate|
|vial_ID_vars	|Integer |	Number of variables in naming convention that are consistent across a time-dependent experiment (ex. genotype, sex)|
|naming_convention |	String |	Experimental conditions in file name|
|path_project| String |	Path to parent folder containing experimental files, configuration_file.cfg and end results.csv file eventually saved here|
|file_suffix |	String |	Suffix of videos being processed|
|convert_to_cm_sec |	Boolean | True if converting output slope to centimeters per second|
