<h1>FreeClimber</h1>

<h3>Overview</h3>

`FreeClimber` is a Python 3-based, background-subtracting particle detection algorithm that performs a local linear regression to quantify the vertical velocity of points moving in a common direction. 

<img src="https://github.com/adamspierer/FreeClimber/blob/master/z/0-Tutorial_climbing.gif" width="600" height="200" align="center">

In lay terms:

- **Background-subtracting**: Removes background pixels that do not change in a set range of frames over a video.

- **Particle detection**: Identifies the x,y-coordinates of each spot/point/marker by frame/time. 

- **Local linear regression**: Finds the most linear segment of a position vs. time (velocity) curve by testing for the greatest regression coefficient over a subset of consecutive frames.

- **Vertical velocity**: Slope of the most linear segment of the velocity curve, which corresponds with the most consistent increase in mean vertical position across n-frames.

This program was designed initially for assessing climbing performance in a *Drosophila* (fruit fly) [rapid iterative negative geotaxis (RING) assay](https://www.sciencedirect.com/science/article/pii/S0531556505000343?casa_token=E8QE2aYrEwoAAAAA:MNa-Wc8BeOXMvmlNuj-b4tH2cMQFuI1ZfUt8qZm0IRY8Qe88xOvw0em07UpwkNqh0QBIPbNZikY). However, this program employs several functions that may be of use beyond the initial design and can be accessed from the source code. This program includes a Graphical User Interface (GUI) for optimizing parameter configurations, and a command line interface for batch processing. This platform has several advantages over over methods and circumvents systemic biases associated with manual methods that are traditionally used to quantify climbing performance in flies.

<h3>Requirements</h3>

    - argparse   [1.1   ]
    - ffmpeg     [1.4   ]
    - matplotlib [3.1.3 ]
    - numpy      [1.18.1]
    - pandas     [1.0.0 ]
    - pip        [20.0.2]
    - scipy      [1.4.1 ]
    - trackpy    [0.4.2 ]
    - wxPython   [4.0.4 ]

NOTE: This platform was developed using Python 3.7.6, and also tested in a Python 3.6 environment. Earlier Python 3 versions will likely work as well but are untested.

<h3>Installing</h3>

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

**Using `conda`**:

	conda install -c adamspierer freeclimber

**Using PyPi**: 

	pip install FreeClimber

**Download the script files** (can be done with `git clone` if user is familiar with `git` or by directly downloading the `.py` files into a single folder.

**Cloning the git repository**:

	cd <folder of interest>
	git clone https://github.com/adamspierer/FreeClimber.git
	
NOTE: As of now, the platform itself is <u>not</u> a module and these steps merely download the dependencies. The script files must be directly referenced when running the program. See our [tutorial](<link to tutorial file>) for usage instructions.


<h3>Test files</h3>

There are three different sets of video files in the `test` folder to demonstrate varying conditions.

| `clean_background` | `diff_number_flies` | `noisy_background` | 
| --- |---| ---|
| <img src="https://github.com/adamspierer/FreeClimber/blob/master/z/clean_setup.jpg" width="366" height="200">       | <img src="https://github.com/adamspierer/FreeClimber/blob/master/z/diff_number_flies.jpg" width="366" height="200"> |  <img src="https://github.com/adamspierer/FreeClimber/blob/master/z/noisy_background.png" width="366" height="200">  |

    - clean_background: contains three vials of 9 flies per vial set on visually clean background (./examples/w1118_m_2_1.h264 OR ./examples/other/clean_background_few_flies.h264)
    - diff_number_flies: contains six vials with 5, 10, 19, 23, 15, 33 flies per vial set on visually clean background (./examples_other/clean_background_different_flies.h264)
    - noisy_background: small, noisy background with low resolution. This example is meant to illustrate how the platform can be configured to work with sub-optimal backgrounds, though better videos yield better results. (./examples/noisy_background_few_flies.h264)

We've also included an additional video with a wagon-wheel effect to illustrate how the detector can still perform well, even with horizontal lines moving across the screen. See ./example\_other/clean\_background\_many\_flies\_with\_lines.h264

Video files were recorded in h264 format (Raspberry Pi default) and can be incompatibile with certain media players. We recommend using [VLC media player](https://www.videolan.org/vlc/index.html).

<h3>Usage</h3>

The following is a general overview of the platform usage. For detailed instructions, please see our [tutorial page](https://github.com/adamspierer/FreeClimber/blob/master/TUTORIAL.md).

Make sure the FreeClimber scripts are downloaded and in a folder on your computer. Navigate to the `FreeClimber` directory and type:

	cd <path_to_FreeClimber_scripts>

**To run from the GUI**:

Specify a file path to a video after the `--video_file` flag.

    pythonw ./scripts/FreeClimber_gui.py --video_file ./example/w1118_m_2_1.h264

**To run from the command line**:

Specify a file path to a configuration file after the `--config_file` flag.

    python FreeClimber_main.py --config_file ./example/example.cfg
    
--> A required `<file>.cfg` is needed to run the command line tool. This file is generated by the GUI, or can be modified from the provided example file. The program will find all specified video files of a common type (`file_suffix`) that are nested within the specified parent folder (`project_path`). 

There are other options for running the command line version that can be accessed using the `-h` flag. In particular, there are options to `--process_all`, `--process_undone`, and `--process_custom`. As the names imply, you can process all, those unprocessed, or a custom list of flies. See the [Tutorial page]('https://github.com/adamspierer/FreeClimber/blob/master/TUTORIAL.md') for more information.

<h3>Code Structure/Overview</h3>

`FreeClimber_gui.py` - GUI wrapper for the detector.py script.

`FreeClimber_main.py` - Command line interface wrapper for the detector.py script

`detector.py` - Contains the detector object, which is important for parsing the video file into a multi-dimensional numpy array (ndarray) and all the functions needed to get the data out.

`example.cfg` - A configuration file generated by the GUI or modified by the user from the provided example. It contains the detection parameters needed to run the program.

`video_file.suffix` - A video file (not image stack) to read in. Unsupported file formats can be converted using [FFmpeg](https://www.ffmpeg.org/), a well-documented multimedia editing platform.

`gather_files.py` - Outputs a list of all files with a given suffix that can be used for customizing the files a user wants to process.

`custom.prc` - Output from `gather_files.py`, contains file paths that can be used to customize the files FreeClimber processes.

We encourage you to to visit our [Tutorial page]('https://github.com/adamspierer/FreeClimber/blob/master/TUTORIAL.md') for a more thorough walk-through, description, and various caveats.

<h3>Deployment</h3>

This software has only been tested on a Mac OS X (Sierra 10.12.6) but is likely not limited to this OS.


<h3>Contributing</h3>

Contributors can fork from the repository and submit a pull request when modifications are ready. Please document the changes you made and any pertinent information that will help in our review of the changes.


<h3>Release History</h3>

We plan to release maintenance updates as needed, though we are unlikely to modify the platform's main functionality.


<h3>Citing this work</h3>
The manuscript associated with this platform is in the end stages of revision or in the review process. Please contact the authors directly for how to cite this work.

<h3>License</h3>

This work is licensed under the MIT license.

<h3>Authors</h3>

Written by [Adam Spierer](https://github.com/adamspierer) and [Lei Zhuo](https://github.com/ctzhu/) with special thanks to Brown University's [Computational Biology Core](https://github.com/compbiocore/) for assistance with code review.