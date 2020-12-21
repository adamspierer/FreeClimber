<h1>FreeClimber</h1>

<h3>Overview</h3>

`FreeClimber` is a Python-based, background-subtracting particle detection algorithm that performs a local linear regression to quantify the vertical velocity of points moving in a common direction.

<img src="https://github.com/adamspierer/FreeClimber/blob/master/z/0-Tutorial_climbing.gif" width="600" height="200" align="center">

In lay terms:

- **Background-subtracting**: Removes a video's background pixels to improve particle detection.

- **Particle detection**: Identifies the x,y,time-coordinates of each spot/marker. 

- **Local linear regression**: Finds the most linear segment of a position vs. time (velocity) curve, via linear regression over each subset of n-consecutive frames.

- **Vertical velocity**: Slope of the most linear segment of the velocity curve, which corresponds with the most consistent increase in mean vertical position across n-frames.

This program was designed initially for assessing climbing performance in a *Drosophila* (fruit fly) [rapid iterative negative geotaxis (RING) assay](https://www.sciencedirect.com/science/article/pii/S0531556505000343?casa_token=E8QE2aYrEwoAAAAA:MNa-Wc8BeOXMvmlNuj-b4tH2cMQFuI1ZfUt8qZm0IRY8Qe88xOvw0em07UpwkNqh0QBIPbNZikY). However, this program employs several functions that may be of use beyond the initial design and can be accessed from the source code or loaded as a module. This program includes a Graphical User Interface (GUI) for optimizing parameter configurations, and a command line interface for batch processing. This platform has several advantages over over methods and circumvents systemic biases associated with manual methods that are traditionally used to quantify climbing performance in flies.

<h3>Requirements</h3>

General programs:
	- FFmpeg        [4.3.1 ]

Python modules:
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

**Installation using PyPi**

	pip install FreeClimber

<!---
**Using conda (still in testing)**:

    conda install -c adamspierer freeclimber
-->

**Download the script files** (can be done with `git clone` if user is familiar with `git` or by directly downloading the `.py` files into a single folder)

**Cloning the git repository**:

	cd <folder of interest>
	git clone https://github.com/adamspierer/FreeClimber.git
	
NOTE: As of now, the platform itself is <u>not</u> callable as a module and these steps merely download the dependencies. The script files must be directly referenced when running the program. See our [tutorial](https://github.com/adamspierer/FreeClimber/blob/master/TUTORIAL.md) for usage instructions.

<h3>Test files</h3>

The `example` folder contains the video file used in the Tutorial. These files contain the video and resulting data and plot files are also included.

**Inputs** (by file suffix):

In ./example folder (main test file):

- **.h264** - Video file (default Raspberry Pi output, see note below).
- **.mov** - Video file (more standard video file type).
- **.cfg** - Configuration file (generated from GUI but can be modified like a text document).

In ./example\_other folder (additional test files):

- **/ex_1/climbing_1.h264** - Clean background, rig from example folder, with different numbers of flies per vial.
- **/ex_2/climbing_2.h264** - Clean background, rig from example folder, with similar numbers of flies per vial, plus wagon-wheel banding effect.
- **/ex_3/climbing_3.h264** - Clean background, homemade rig 1, with similar numbers of flies per vial. Same video as gif in overview.
- **/ex_4/climbing_4.h264** - Heterogeneous background, homemade rig 2/test tube holder, with similar numbers of flies.
- **/ex_5/climbing_5.h264** - Heterogeneous background, homemade rig 3, borderline unusable quality.

Outputs (by file suffix):

- **.raw.csv** - Data - All spots before filtering. Can be used as input with trackpy for tracking spots.
- **.filtered.csv** - Data - Spots after filtering, binning, and labeling.
- **.slopes.csv** - Data - Each vial's slope and best local linear regression score.
- **.ROI.png** - Plot - First frame of video with a box drawn around the ROI, white dividers corresponding with vial bins, and a cyan box drawn around the ROI post-edge filtering.
- **.processed.csv** - Plot - Three subplots corresponding with the cropped and grayscaled frame (variable `frame_0`), null background image, and background subtracted image.
- **.spot_check.png** - Plot - Three sets of two subplots. Sets correspond with three filtering parameters: eccentricity (ecc, circularity (0 circular <--> 1 not circular), spot mass, and spot signal. The two subplots for each correspond with a histogram colored according to the spot metric value on the corresponding scatterplot for the x,y location of each spot in the video (all overlayed on the first frame of the video).
- **.diagnostic.png** - Plot - Overlay scatterplots for the spots identified in the first and last frame of the local linear regression for all spots in the video, the mean-vertical position vs. time plot for each vial (darker section corresponds with most linear), and an overlay scatterplot for all points identified throughout the video.

Other outputs:

- **results.csv** - Merging of all video .slopes.csv files into a single file.
- **log/completed.log** - File paths for videos fully processed.
- **log/skipped.log** - File paths for videos skipped during processing.


Video files were recorded in .h264 format (Raspberry Pi default) and can be incompatibile with certain media players. We recommend using [VLC media player](https://www.videolan.org/vlc/index.html), though provide the .mp4 format as well. Nearly all common video file types should be compatible with FreeClimber, though ones that are not can be converted using [FFmpeg](https://ffmpeg.org/).

<h3>Usage</h3>

The following is a general overview of the platform usage. For detailed instructions, please see our [tutorial page](https://github.com/adamspierer/FreeClimber/blob/master/TUTORIAL.md).

Make sure the FreeClimber scripts are downloaded and in a folder on your computer. Navigate to the `FreeClimber` directory and type:

	cd <path_to_FreeClimber>

**To run from the GUI**:

Specify a file path to a video after the `--video_file` flag.

    pythonw ./scripts/FreeClimber_gui.py --video_file ./example/<video_file.suffix>

- Example video_file.suffix: `w1118_m_2_1.h264`

**To run from the command line**:

Specify a file path to a configuration file after the `--config_file` flag.

    python FreeClimber_main.py --config_file ./example/example.cfg
    
- A required `<file>.cfg` is needed to run the command line tool. This file is generated by the GUI, or can be modified from the provided `example.cfg` file. The program will find all specified video files of a common type (`file_suffix`) that are nested within the specified parent folder (`project_path`), rather than where the configuration file is held. 

- There are other options for running the command line version that can be accessed using the `-h` flag. Notably, there are options for processing `--process_all`, `--process_undone`, and `--process_custom`. As the names imply, you can process all, those unprocessed, or a custom list of flies. See the [Tutorial page]('https://github.com/adamspierer/FreeClimber/blob/master/TUTORIAL.md') for more information.

<h3>Code Structure/Overview</h3>

`FreeClimber_gui.py` - GUI wrapper for the detector.py script.

`FreeClimber_main.py` - Command line interface wrapper for the detector.py script.

`detector.py` - Contains the detector object, which is important for parsing the video file into a multi-dimensional numpy array (ndarray) and all the functions needed to get the data out.

`example.cfg` - A configuration file generated by the GUI or modified by the user from the provided example. It contains the detection parameters needed to run the program.

`gather_files.py` - Outputs a list of all files with a given suffix that can be used for customizing the files a user wants to process.

`custom.prc` - Output from `gather_files.py`, contains file paths that can be used to customize the files FreeClimber processes. File paths from the `log/skipped.log` file can also be copied and pasted into a similar file.

We encourage you to to visit our [Tutorial page]('https://github.com/adamspierer/FreeClimber/blob/master/TUTORIAL.md') for a more thorough walk-through, description, and various caveats.

<h3>Version releases</h3>
0.4.0 - Fixed error with frame cropping and added additional variable compatibility checks.

0.3.2 - Release version for Journal of Experimental Biology

<h3>Deployment</h3>

This software has only been tested on a Mac OS X (Sierra 10.12.6) but is likely not limited to this OS.

<h3>Contributing</h3>

Contributors can fork from the repository and submit a pull request when modifications are ready. Please document the changes you made and any pertinent information that will help in our review of the changes.

<h3>Release History</h3>

We plan to release maintenance updates as needed, though we are unlikely to modify the platform's main functionality.

<h3>Citing this work</h3>
If you use FreeClimber in your work, please cite:
> A.N. Spierer, D. Yoon, CT. Zhu, and DM. Rand. (2020) FreeClimber: Automated quantification of climbing performance in <i>Drosophila</i>. <i>Journal of Experimental Biology</i>. DOI: 10.1242/jeb.229377

<h3>License</h3>

This work is licensed under the MIT license.

<h3>Authors</h3>

Written by [Adam Spierer](https://github.com/adamspierer) and [Lei Zhuo](https://github.com/ctzhu/) with special thanks to Brown University's [Computational Biology Core](https://github.com/compbiocore/) for assistance with code review.