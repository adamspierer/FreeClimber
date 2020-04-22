<h1>FreeClimber -- Python 3</h1>

<h3>Overview</h3>

`FreeClimber` is a Python-based, particle detection algorithm used for quantifying the climbing speed for a group of flies. This platform has several advantages over current methods. First, its 's advantages include its use of a background subtraction step to clean otherwise noisy backgrounds, as well as a sliding window to calculate the most linear portion of a mean y-position vs. time curve to find the most accurate representation of a cohort's climbing velocity. This gif is a good demonstration of that, where the most linear increase in position over a 2-second window is blue--with an accompanying dashed line of best fit--while the remaining points are red. 

<img src="https://github.com/adamspierer/vial_detector/blob/python3/supplemental/tutorial_0.gif" width="600" height="200" align="center">

While this platform is most useful for quantifying a negative geotaxis (climbing) assay with Drosophila, it can be used for any video with a static background, and where the particles are of approximately uniform size and traveling from the bottom of the frame to the top.
 


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

NOTE: This platform was developed using Python 3.7.6, and while untested, earlier versions of Python 3 will likely work as well.

<h3>Installing</h3>
<ol><li>We recommend running this package in an Anaconda [download](https://docs.anaconda.com/anaconda/install/) virtual environment running Python 3 [how to](https://uoa-eresearch.github.io/eresearch-cookbook/recipe/2014/11/20/conda/).</li>

<li>The aforementioned requirements can be downloaded using conda:
```conda install -c adamspierer freeclimber```
or PyPi:
```pip install FreeClimber```</li>

<ul><li>NOTE: Some systems may have issues with dependencies in a Python 3.7 virtual environment, so we recommend a Python 3.6 environment if this is the case.</li></ul>

<li>Clone the git repository, or simply download the files with the `.py` suffix to a single folder.</li></ol>

<h3>Test files</h3>

There are three different sets of video files in the `test` folder to demonstrate varying conditions.

| clean_background | diff_number_flies | noisy_background | 
| --- |---| ---|
| <img src="https://github.com/adamspierer/vial_detector/blob/python3/supplemental/example_images/clean_setup.jpg" width="366" height="200">       | <img src="https://github.com/adamspierer/vial_detector/blob/python3/supplemental/example_images/diff_number_flies.jpg" width="366" height="200"> |  <img src="https://github.com/adamspierer/vial_detector/blob/python3/supplemental/example_images/noisy_background.png" width="366" height="200">  |

    - clean_background: contains three vials of 9 flies per vial set on visually clean background
    - diff_number_flies: contains six vials with 5, 10, 19, 23, 15, 33 flies per vial set on visually clean background
    - noisy_background: small, noisy background with low resolution. This example is meant to illustrate how the platform can be configured to work with sub-optimal backgrounds, though better videos yield better results.

<h3>Usage</h3>

Navigate to the `FreeClimber` directory:

```cd <path_to_FreeClimber_scripts>```

To run the GUI, navigate to working folder and type:

    ```pythonw FreeClimber_gui.py <optional: path_to_video_file>```

--> Optional argument will take a video file path and automatically load it. Leaving this blank will cause a dialog box to open for the user to search for a video to load.

To run from the command line, navigate to the working folder and type:

    ```python FreeClimber_main.py <path_to_configuration_file.cfg> <optional: False>```

--> A required `<file>.cfg` is needed to run the command line tool. This file is generated by the GUI, or can be modified from the provided example file. The program will process all files with the specified file suffix (ex. h264, mov, avi) in the `path_project` folder.

We provide a more complete tutorial on how to run this program: https://github.com/adamspierer/vial_detector/blob/python3/supplemental/Tutorial.md

<h3>Code Structure/Overview</h3>

`FreeClimber_gui.py` - Sets up the GUI object

`DectFrame.py` - Organizes and formats the GUI

`detector.py` - Contains the *detector* object, which is important for parsing the video file into a multi-dimensional numpy array (ndarray) and all the functions needed to get the data out.

`FreeClimber_main.py` - This is the main script; while takes in a configuration file with the necessary detection parameters and processes videos.

`example.cfg` - A configuration file generated by the GUI or modified by the user. It contains the detection parameters needed to run the program.

`video_file.suffix` - A video file to read in. Currently works with most common movie formats (pre-approved suffixes include: .flv, .h264, .mov, .mp4, .wmv, and more). There is a known issue with .avi, so convert to a different format using:

```ffmpeg -i <your_file>.avi -c copy <your_file>.<pre-approved suffix>```

We encourage you to to visit our [Tutorial page]('https://github.com/adamspierer/vial_detector/blob/python3/supplemental/Tutorial.md') to explore this program further.


<h3>Deployment</h3>

This software has only been tested on a Mac OS X (Sierra 10.12.6).


<h3>Contributing</h3>

This will be the only release of this program.


<h3>Release History</h3>

This will be the only release of this program.


<h3>Citing this work</h3>
<Citation for this paper>

<h3>Authors</h3>

Written by Adam Spierer and Lei Zhuo with special thanks to Brown University's Computational Biology Core for assistance with code review.