<h1>Free Climber</h1>

<h3>Overview</h3>

This tutorial will walk you through how to use the `free_climber` platform. As of now, you must specify the path to the `free_climber` script you want to run, rather than running from a module. The goal is to make this all a module moving forward, but for this stage in development here's what I've got...

The `free_climber` platform is designed to take a video of flies climbing in vials, and quantify the group's velocity. This gif demonstrates the principle underlying how velocity is quantified:

<img src="https://github.com/adamspierer/free_climber/blob/python3/supplemental/tutorial_0.gif" width="500" height="167">

Blue '+' represent candidate spots, while colored circles represent 'true' spots that passed a cutoff threshold. Blue circles represent the portion of the video where the mean cohort position is rising most consistently (seen in both the video of flies climbing and in the position vs. time plot). A black dotted line signifies the line of best fit through this portion of the curve, while a grey jagged line represents the trajectory of points throughout the video. As flies reach the top of the vial, they become more challenging for the detector to find so the curve begins to plateau and becomes less linear. When this happens, the blue circles become red and the final climbing metric is no longer calculated.

This platform enables users to use an interactive, graphical user interface (GUI) to optimize detection parameters and specify experimental details, as well as an option for batch processing many videos using the command line. This tutorial will walk through both approaches.

<h2>Installation and requirements</h2>

To install, make sure `conda` is [installed on your computer](https://docs.anaconda.com/anaconda/install/).

Then, type: `conda install -c adamspierer free_climber` into the command line to install the package requirements.

Create a [virtual environment](https://uoa-eresearch.github.io/eresearch-cookbook/recipe/2014/11/20/conda/) for Python3 using conda, since installing some modules directly to the base operating system can cause errors.

Platform was written in Python 3.7.6 and requires the following modules:

    - matplotlib [3.1.3 ]
    - numpy      [1.18.1]
    - pandas     [1.0.0 ]
    - pip        [20.0.2]
    - python     [3.7.6 ]
    - scipy      [1.4.1 ]
    - trackpy    [0.4.2 ]
    - wxPython   [4.0.4 ]

Also requires:
    - ffmpeg (only tested with v.4.0, but should run with others)

Alternatively, the Github repository can be cloned and dependencies manually installed.

<h2>Running the interactive, Graphical User Interface</h2>

The graphical user interface (GUI) is useful for optimizing the parameters one might use in analyzing a video. Several fields are available to modify and buttons enable quick testing.

<h3>Step 1 - Specify a video</h3>

Assuming this is the first walkthrough, we should run through the GUI to optimize the detection parameters. To begin, run `pythonw free_climber_gui.py`. A dialog box will appear and you should select a video to start with. If you picked the wrong file, you can always reselect using the `Browse...` button in the upper left.

When the video loads, the backend will scrape some of the metadata from the video and populate some of the fields. One of these fields is the `Frames / sec` field.

**Frame rate** - The default value is 25, so if the frame rate is known but different, then manually override it or leave it and keep it consistent across videos.

<h3>Step 2 - Defining video parameters</h3>

Selecting-and-dragging the cursor across the video will create a red rectangle will appear. This rectangle can be redrawn several times, making it useful for measuring known distances to calibrate the final slopes into pixels per centimeter and for drawing a Region of Interest (ROI) around the vials.

**Pixels per centimeter** - To calculate the number of pixels per centimeter, draw a horizontal or vertical line between two points of known distance (e.g. the width of a fly vial -- narrow glass = 2.45 cm across). Record the corresponding value in the `Width` or `Height` and resample if there are other markers in frame. Take the average of these measurements and divide by their known distance in cm. Enter this new value into the `Pixels / cm` field. 

For example, if three vials were measured at 110, 114, 111 -- then the average would be 111.7, which is then divided by 2.45 to equal 45.6 pixels/cm. If you want the final slope values output from the platform to be in cm / second, then click the `Convert to cm?` box directly below the field you populated. Otherwise, resulting slope values will be in pixels / second.

<img src="https://github.com/adamspierer/free_climber/blob/python3/supplemental/2-Tutorial_measure_vials.png" width="500" height="400">        <img src="https://github.com/adamspierer/free_climber/blob/python3/supplemental/pixel_to_cm_calculation.png" width="223" height="200">

**Draw a Region of Interest (ROI)** - Move the cursor over to the upper left-most point, press and hold as you draw a rectangle over the Region Of Interest (ROI) you want to analyze. It is important that the **rectangle be drawn from the upper left to the lower right**--negative values for the height and/or width are not supported. Some may want to draw a larger rectangle if there is variation in the vial rack's position between videos, while others may want a rectangle that is tighter to the experimental rig if there is little positional variation between videos. The larger the box, the more pixels to process, the longer it will take to process a video. 

<img src="https://github.com/adamspierer/free_climber/blob/python3/supplemental/3-Tutorial_draw_box.png" width="500" height="400">

**Fixed ROI size?** - This box allows the user to specify a custom `width` and `height` for the rectangle, and then click anywhere on the image and have that box move around.


<h3>Step 3 - Specify spot parameters</h3>

These parameters are important for identifying candidate spots, as well as setting a threshold for determining whether a spot is "True" or not.

**Diameter** - Approximate diameter for estimating a spot. Only odd numbers, rounding up when in doubt.

**MaxDiameter** - The maximum distance across to consider a spot.

**MinMass** - The minimum integrated brightness.

**Threshold** - The threshold value for calling a spot "True" or "False." True/False spots are visualized using the `Test Parameters` button, more below.
For more information, see the [TrackPy documentation](http://soft-matter.github.io/trackpy/v0.4.2/generated/trackpy.locate.html#trackpy.locate). Start low for the first pass at optimizing the parameters--the `Test parameters` button will output a plot that helps better gauge an appropriate value later on.

<h3>Step 4 - Additional parameters</h3>

**Background Frames** - The first through last frames from which to calculate a background image. The default is to calculate it from the entire video, but a subset can be taken to insignificantly reduce processing time. The median pixel intensity across all coordinates is calculated (middle image), and when that is subtracted from an individual frame (left image), the program is left with just the pixels that change color (signifying movement; right image). Here is an example of what this looks like to the computer:

<img src="https://github.com/adamspierer/free_climber/blob/python3/supplemental/background_subtraction/4_compilation.png" width="600" height="200">

**Check Frames** - Enter two separate frames to display the candidate and True points when running the `Test Parameters` button. The defaults are the first frame and the one calculated to be the frame at 2 seconds.

**Number of vials** - The number of vials to bin flies into within a video. Vials are evenly spaced bins across the minimum and maximum x-coordinates in a video. This is useful for following unique vials across time or treatments. That said, it is important to have flies in the edge vials since the minimum and maximum x-coordinates in a video determine how the spots are binned across vials. Alternatively, videos can be analyzed as a single vial by entering 1.

**Window size** - Number of frames for sliding window when calculating the local linear regression. The best practice is to input the number of frames corresponding with the most linear portion of the sigmoidal, mean y-position vs. time plot and use a few frames short of the length of that linear portion. We found flies can climb the vial in 1 to 2 seconds. The longer the length, the better the estimate. But if the length is too long it may be a poor estimate. 

<h3>Step 5 - Naming parameters</h3>

**Naming pattern** - All input files should be named with experimental details separated by an underscore. For instance, if an experiment was testing genotype and sex performance across days, with several replicates, then the naming convention might be `genotype_sex_day_rep`. The order is somewhat important; experimental details consistent during the experiment (e.g. genotype and sex) should be in the beginning, while those that change (e.g. time and replicate) should be at the end. The order of variables within those groupings does not matter.

**Video vars** - The number of variables that are consistent during the experiment (e.g. 2 for genotype and sex). See entry above.

**Project path** - The parent folder to search through for all video files of a given type. The `.cfg` file and `results.csv` file containing all slopes from all videos at the end will eventually be saved to this folder.


<h3>Buttons</h3>

**Test parameters** - This button allows the user to optimizing the spot parameters and background frames. All "false" spots have a blue '+' over them, while "true" spots have a red circle over them. "True" vs "false" spots are labeled such if the 'signal' for that spot exceeds the user-defined `Threshold`. A histogram plot will assist in determining an appropriate `Threshold`. Just look for the region just before the large hump on the far right. After resulting plots are generated, one can reload the same video with the `Reload video`, or press any of the other buttons.

<img src="https://github.com/adamspierer/free_climber/blob/python3/supplemental/4-Tutorial_test_parameters.png" width="500" height="400">

**Reload video** - Reload a video to refine detection parameters once the `Test parameters` button has been pressed.

**Save parameters** - Generate a `.cfg` file from the displayed parameters. This file can be edited manually to use the same analysis parameters on a different set of files.

**Run analysis** - Perform full analysis on the selected video.

<img src="https://github.com/adamspierer/free_climber/blob/python3/supplemental/5-Tutorial_analysis.png" width="500" height="400">

<h2>Running the command line version</h2>

To use, run: 

```python ./free_climber_main.py <path_to_configuration_file.cfg>```

If you don't want all `<file_name.slopes.csv>` files concatenated within a project, then run: `python ./free_climber_main.py <path_to_configuration_file.cfg> False`


<h4>Tips</h4>

- We recommend trimming the length of the video down to the length you are interested in analyzing and having a single replicate of the assay per video. To slice videos efficiently using a command line-based tool, we recommend using `FFmpeg`, which is well documented: https://trac.ffmpeg.org/wiki/Seeking.

- For best spot-detection, we also recommend videos have a clean and consistent background. One of the advantages to this platform is the method of background subtraction. As long as the background is static and there is sufficient contrast between it and the flies, the detector should work reasonably well (see the `noisy_background` example). Even so, we found backlighting the flies with an LED-tracing board created strong contrast and made for better data quality. 

- If using an LED-tracking board, one issue that can arise is a 'wagon-wheel effect' where a "rolling shutter" paired with a repetitive/cyclical phenomenon creates an optical illusion. In this case, the LEDs will create a horizontal banding effect that will affect the ability of the detector to subtract the background and negatively affect the data quality. Adjust the dimmer settings until you find a brightness that works. **GIF TO COME**

- When optimizing spot detection parameters, first focus on the `Diameter`, `MaxDiameter`, and `MinMass`. These fields determine whether a group of pixels is a candidate spot so find the highest values for these parameters that cover all obvious spots, since lower values can increase the time required for the program to process a video. The `Threshold` field should be optimized last to determine which spots pass a signal threshold. A hint here is to find a `Threshold` value that covers most of the points in `Check frames` when after `Test parameters`. Open the <file-name.raw.csv> file and find the frames' spot ST_dsty values to find inform the `Threshold` field. 

- Image coordinates are numbered from upper left to lower right, so <file_name.raw.csv> and <file_name.filter.csv> will have y-coordinates inverted for plotting purposes. This is important for those looking to use the outputs for further visualization.

- The `<file_name.raw.csv>` output is also really useful for passing to `TrackPy` for linking spots together into tracks. Read more on TrackPy particle __tracking__ [API](http://soft-matter.github.io/trackpy/v0.4.2/api.html). Briefly, the file can be loaded in and saved as `f`, which will seamlessly integrate into the `Step 3: Link features into particle trajectories` in the [TrackPy Walkthrough](http://soft-matter.github.io/trackpy/v0.4.2/tutorial/walkthrough.html).

<h4>Table of variables</h4>
|Variable name |	Data type |	Explanation|
| --- | --- | --- |
|x |	Integer |	Leftmost pixel of ROI for analysis|
|y |	Integer	| Topmost pixel of ROI for analysis|
|w | Integer |	Width of ROI for analysis|
|h | Integer |	Height of ROI for analysis|
|frame_0 |	Integer	| First frame to display in Test Parameters|
|frame_n |	Integer	| Last frame to display in Test Parameters|
|blank_0 |	Integer	| First frame of range to subtract background|
|blank_n |	Integer	| Last frame of range to subtract background|
|threshold |	Integer	| Threshold for filtering against points|
|maxsize | Integer| The maximum distance across to consider a spot|
|minmass	|Integer	| The minimum integrated brightness|
|diameter	|Integer	| Estimated diameter of spot in pixels|
|vials	| Integer |	Number of vials in video|
|window |	Integer	| Number of frames for sliding window|
|pixel_to_cm|	Integer	| Conversion factor for pixels to centimeters|
|frame_rate	|Integer |	Video frame rate|
|vvars	|Integer |	Number of variables in naming convention that are consistent across a time-dependent experiment (ex. genotype, sex)|
|convert_to_cm |	Boolean | True if converting output slope to centimeters per second|
|naming_convention |	String |	Experimental conditions in file name|
|path_project| String |	Path to parent folder containing experimental files, configuration_file.cfg and end results.csv file eventually saved here|
|file_suffix |	String |	Suffix of videos being processed|