#!/usr/bin/env python
# -*- coding: utf-8 -*- 

## File name : detector.py
## Created by: Adam N. Spierer
## Date      : December 2020
## Purpose   : Script contains main functions used in FreeClimber package, as well as added functionality

version = '0.4.0'
publication = False

import os
import sys
import time
import ffmpeg

import numpy as np
import pandas as pd
import trackpy as tp
import subprocess as sp
from scipy.stats import linregress
from scipy.signal import find_peaks,peak_prominences

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D

## Issue with 'SettingWithCopyWarning' in step_3
pd.options.mode.chained_assignment = None  # default='warn'

class detector(object):
    '''Particle detection platform for identifying the group climbing velocity of a 
    group of flies (or particles) in a Drosophila negative geotaxis (climbing) assay.
    This platform is designed to take videos of varying background homogeneity and 
    applying a background subtraction step before extracting the x,y,time coordinates of
    spots/flies. Climbing velocities are calculated as the most linear portion of a user-defined
    subset of frames by vial (vertical divisions of evenly spaced bins from the min/max
    X-range.
    '''
    def __init__(self, video_file, config_file = None, gui = False, variables = None, debug = False, **kwargs):
        '''Initializing detector object
        ----
        Inputs:
          video_file (str): Path to video file to process
          config_file (str): Path to configuration (.cfg) file
          gui (bool): GUI-specific functions
          variables: None if importing from a file, or list if doing so manually
          debug (bool): Prints out each function as it runs.
          **kwargs: Keyword arguments that are unspecified but can be passed to various plot functions
        ----
        Returns:
          None -- variables are saved to the detector object
        '''
        self.debug = debug
        if self.debug: print('detector.__init__')
        
        self.config_file = config_file
        self.video_file = self.check_video(video_file)
        
        ## Load variables
        if gui:
            self.load_for_gui(variables = variables)
        elif  not gui:
            self.load_for_main(config_file = config_file)
        self.check_variable_formats()

        ## Setting a color map
        self.vial_color_map = cm.jet
        
        ## Create a conversion factor
        if self.convert_to_cm_sec: self.conversion_factor = self.pixel_to_cm / self.frame_rate
        else: self.conversion_factor = 1

        print('')
        self.specify_paths_details(video_file)
        self.image_stack = self.video_to_array(video_file,loglevel='panic')
        return

    ## Loading functions    
    def load_for_gui(self,variables):
        '''Loads experimental and detections variables to the detector object for GUI
        ----
        Inputs:
          variables (list): List of variables found in the configuration file, order is irrelevant
        ----
        Returns:
          None -- variables are saved to the detector object'''    
        if self.debug: print('detector.load_for_gui')
        self.config = None
        self.path_project = None
        
        ## Exit program if no variables are passed
        if variables == None:
            print('\n\nExiting program. No variable list loaded')
            raise SystemExit
        
        ## Pass imported variables to the detector object 
        if self.debug: print('detector.load_for_gui: --------variables--------')
        for item in variables:
            if self.debug: print('detector.load_for_gui:',item)
            if ~item.startswith(('\s','\t','\n')):
                try: exec('self.'+item)
                except: print('detector.load_for_gui: !! Could not import ( %s )' % item) 
        return     
        
    def load_for_main(self, config_file = None):
        '''Loads experimental and detections variables to the detector object for command 
        line interface.
        ----
        Inputs:
          config_file (str): Path to configuration (.cfg) file
        ----
        Returns:
          None -- variables are saved to the detector object''' 
        if self.debug: print('detector.load_for_main')
        
        ## Read in lines from configuration (.cfg) file
        try:
            with open(config_file,'r') as f:
                variables = f.readlines()
            f.close()
            
            ## Filter, format, and import variables to detector object
            if self.debug: print('detector.load_for_main:  --------variables--------')
            variables = [item.rstrip() for item in variables if not item.startswith(('#','\s','\t','\n'))]

            for item in variables:
                if self.debug: print('detector.load_for_main:',item)
                if ~item.startswith(('\s','\t','\n')):
                    try: exec('self.'+item)
                    except: print('detector.load_for_main: !! Could not import ( %s )' % item) 
            return

        ## Exit program if issue with the configuration file
        except: 
            print('\n\nExiting program. Could not read in file.cfg, but path and suffix are good--likely a formatting issue')
            raise SystemExit
        return

    ## Checking video path is valid
    def check_video(self,video_file=None):
        '''Checking video path is valid, exiting if not.
        Input:
          video_file (str): Video file path
        ----
        Returns:
          video_file (str): Video file path'''
        if self.debug: print('detector.check_video...', end='')
        
        ## Check if file path is valid, exit if not
        if os.path.isfile(video_file):
            return video_file
        else:
            print('\n\nExiting program. Invalid path to video file: ',video_file)
            raise SystemExit

    ## Specifying variables
    def specify_paths_details(self,video_file):
        '''Takes in the video file and other imported variables and parses them as needed.
        ----
        Inputs:
          video_file (str): Video file path
        ----
        Returns:
          None -- Passes parsed variables back to detector object
        '''
        if self.debug: print('detector.specify_paths_details')
        
        ## Set file and path names
        folder,name = os.path.split(video_file)
        self.name = name[:-5]
        self.name_nosuffix = '.'.join(video_file.split('.')[:-1])
        
        ## Defining final file names and destinations
        file_names = ['data','filtered','diagnostic','slope']
        file_suffixes = ['.raw.csv','.filtered.csv','.diagnostic.png','.slopes.csv']
        for item,jtem in zip(file_names,file_suffixes):
            var_name = 'self.path_'+item
            file_path = ''.join([self.name_nosuffix,jtem])
            if self.debug: print('detector.specify_paths_details: ' + var_name+"='"+file_path+"'")
            exec(var_name+"='"+file_path+"'")

        ## Project folder specific paths
        if self.path_project == None: self.path_project = os.path.join(folder,self.name + '.cfg')

        ## For future release
#         self.path_review_diagnostic = self.path_project + '/review_at_R_%s/review.log' %self.review_R
        
        ## Extracting file details and naming individual vials
        self.file_details = dict(zip(self.naming_convention.split('_'),self.name.split('_')))
        self.experiment_details = self.name.split('_')
        self.experiment_details[-1] = '.'.join(self.experiment_details[-1].split('.')[:-1])
        self.vial_ID = self.experiment_details[:self.vial_id_vars]
        
        ## Creating a list of colors for plotting
        self.color_list = [self.vial_color_map(i) for i in np.linspace(0,1,self.vials)]
        return

    ## Checking to make sure variables are entered properly...still more to include
    def check_variable_formats(self):
        '''Checks to make sure at least some of the variables input formatted properly'''
        if self.debug: print('detector.check_variable_formats')
        
        ## Vial number
        if self.vials < 1: 
            print('!! Issue with vials: now = 1')
            self.vials = 1
        
        ## Spot diameter must be odd number
        if ~self.diameter%2:
            print('!! Issue with diameter: was %s, now %s' %(self.diameter,self.diameter+1))
            self.diameter += 1
            
        ## Frame rate must be greater than 0
        if self.frame_rate <= 0:
            print('!! Issue with frame_rate: was %s, now 1' %(self.frame_rate))
            self.frame_rate = 1
        
        ## Background blank cannot be greater than the frame 
        if self.blank_0 < self.crop_0:
            self.blank_0 = self.crop_0
        
        if self.blank_n > self.crop_n:
            self.blank_n = self.crop_n

        ## Window size vs. frames to test
        if (self.crop_n - self.crop_0) < self.window:
#             print('!! Issue with window size (%s) being greater than frames (%s). Window size set to 80 percent of desired frames (%s)' %(self.window,self.crop_n - self.crop_0,(.8 * (self.crop_n - self.crop_0))))
            self.window = (self.crop_n - self.crop_0) * 0.8
		
		## blank vs. crop frames
        if self.blank_0 < self.crop_0:
            print('!! Issue with blank frames vs. crop frames. Setting blank_0 (%s) = crop_n (%s)' % (self.blank_0,self.crop_0))
            self.blank_0 = self.crop_0
        if self.blank_n > self.crop_n:
            print('!! Issue with blank frames vs. crop frames. Setting blank_n (%s) = crop_n (%s)' % (self.blank_n,self.crop_n))
            self.blank_n = self.crop_n
        
        ## Check frame is still valid
        if self.check_frame < self.crop_0:
#             print('!! Issue with check_frame < crop_0 (min. cropped frame). Now, check_frame = crop_0 = %s' %self.check_frame)
            self.check_frame = self.crop_0
        if self.check_frame > self.crop_n:
#             print('!! Issue with check_frame > crop_n (max cropped frame). Now, check_frame = crop_n = %s' %self.check_frame)
            self.check_frame = self.crop_n
        return

    ## Video processing functions
    def video_to_array(self, file, **kwargs):
        '''Converts video into an nd-array using ffmpeg-python module.
        ----
        Inputs:
          file (str): Path to video file
          kwargs: Can be used with the ffmpeg.output() argument.
        ----
        Returns:
          image_stack (nd-array): nd-array of the video
        '''
        if self.debug: print('detector.video_to_array')
        
        ## Extracting video meta-data
        try:
            try:
                probe = ffmpeg.probe(file)
            except:
                print('!! Could not read %s into FreeClimber. Likely due to unacceptable video file or FFmpeg not installed' % file)
                raise SystemExit
            video_info = next(x for x in probe['streams'] if x['codec_type'] == 'video')
            self.width = int(video_info['width'])
            self.height = int(video_info['height'])
        except:
            print('!! Could not read in video file metadata')
        
        ## Converting video to nd-array    
        try:
            out,err = (ffmpeg
                       .input(file)
                       .output('pipe:',format='rawvideo', pix_fmt='rgb24',**kwargs)
                       .run(capture_stdout=True))
            self.n_frames = int(len(out)/self.height/self.width/3)
            image_stack = np.frombuffer(out, np.uint8).reshape([-1, self.height, self.width, 3])
        except:
            print('!! Could not read in video file to an array. Error message (if any):', err)

        return image_stack

    def crop_and_grayscale(self,video_array,
                         x = 0 ,x_max = None,
                         y = 0 ,y_max = None,
                         first_frame = None,
                         last_frame = None,
                         grayscale = True):
        '''Crops imported video array to region of interest and converts it to grayscale
        ----
        Inputs:
          video_array (nd-array): image_stack generated from video_to_array function
          x (int): left-most x-position
          x_max (int): right-most x-position
          y (int): lowest y-position
          y_max (int): highest y-position
          first_frame (int): first frame to include
          last_frame (int): last frame to include
          grayscale (bool): True to convert to gray, False to leave in color. 
                                NOTE: Must be in grayscale for FreeClimber, option is 
                                available for functionality beyond FreeClimber
        ----
        Returns:
          clean_stack (nd-array): Cropped and grayscaled (if indicated) video as nd-array'''
        if self.debug: print('detector.crop_and_grayscale')
    
        ## Conditionals for cropping frames and video length
        if first_frame == None: first_frame = self.crop_0
        if last_frame == None: last_frame = self.crop_n
        if y_max == None: y_max = video_array.shape[2]
        if x_max == None: x_max = video_array.shape[1]
        if self.debug: print('detector.crop_and_grayscale: Cropping from frame %s to %s' % (first_frame,last_frame))
    
        ## Setting only frames and ROI to grayscale
        if grayscale:
            if self.debug: print('detector.crop_and_grayscale: Converting to grayscale & cropping ROI to (%s x %s)' % (x_max-x,y_max-y))
            ch_1 = 0.2989 * video_array[first_frame:last_frame,y : y_max,x : x_max,0]
            ch_2 = 0.5870 * video_array[first_frame:last_frame,y : y_max,x : x_max,1]
            ch_3 = 0.1140 * video_array[first_frame:last_frame,y : y_max,x : x_max,2]
            clean_stack = ch_1.astype(float) + ch_2.astype(float) + ch_3.astype(float)
        
        ## Only cropping, no grayscaling
        else:
            clean_stack = video_array[first_frame:last_frame,y : y_max,x : x_max,:]
            
        if self.debug: print('detector.crop_and_grayscale: Final video array dimensions:',clean_stack.shape)
        return clean_stack
    
    ## Subtract background
    def subtract_background(self,video_array=None):
        '''Generate a null background image and subtract that from each frame
        ----
        Inputs:
          video_array (nd-array): clean_stack generated from crop_and_grayscale
          first_frame (int): First frame to consider for background subtraction
          last_frame (int): Last frame to consider for background subtraction
        ----
        Returns:
          spot_stack (nd-array): Background-subtracted image stack
          background (array): Array containing the pixel intensities for each x,y-coordinate'''
        if self.debug: print('detector.subtract_background')
        
        ## Setting the last frame to the end if None provided
        first_frame = self.blank_0
        last_frame = self.blank_n
                    
        ## Generating a null background image as the median pixel intensity across frames
        background = np.median(video_array[first_frame:last_frame,:,:].astype(float), axis=0).astype(int)
        if self.debug: print('detector.subtract_background: dimensions:', background.shape)
        
        ## Subtracting the null background image from each individual frame
        spot_stack = np.subtract(video_array,background)
        return spot_stack, background   


    ## Plots and views
    def view_ROI(self,image = None, border = True, x0 = 0,x1 = None, y0 = 0,y1 = None,
                 color = 'r', bin_lines = None, **kwargs):
        '''Generates image of the first frame w/a rectangle for the region of interest
        ----
        Inputs:
          image (array): Slice (single frame) of the image_stack nd-array.
          border (bool): True draws a rectangle over the region of interest
          x0 (int): Left-most coordinate
          x1 (int): Right-most coordinate
          y0 (int): Top-most coordinate (should also be the lowest y-value)
          y1 (int): Bottom-most coordinate (should also be the highest y-value)
          color (str): Color corresponding with rectangle color for region of interest
          bin_lines (bool): True if drawing lines between calculated vials
          **kwargs: Arguments for plt.imshow
        ----
        Returns:
          None -- Generates an image saved in step_3()
        '''
        if self.debug: print('detector.view_ROI')

        ## Defaults to first frame of image stack if None specified
        if self.debug: print('detector.view_ROI :: Setting frame')
        if image == None:
            image = self.image_stack[0]
        
        ## Plots the slice of nd-array
        if self.debug: print('detector.view_ROI :: Plotting image')
        plt.imshow(image,cmap=cm.Greys_r, **kwargs)
        
        ## Draws a red rectangle over the region of interest
        if self.debug: print('detector.view_ROI :: Draw ROI')
        if border:
            if x1 == None: x1 = self.image_stack[0].shape[1]
            if y1 == None: y1 = self.image_stack[0].shape[0]
            plt.hlines(y0,x0,x1, color = color, alpha = .7)
            plt.hlines(y1,x0,x1, color = color, alpha = .7)
            plt.vlines(x0,y0,y1, color = color, alpha = .7)
            plt.vlines(x1,y0,y1, color = color, alpha = .7, label='ROI')

        ## Draws box to denote where outlier trim lines are
        if self.debug: print('detector.view_ROI :: Trim outlier lines (if selected)')
        if self.trim_outliers:
            lc,rc,tc,bc = self.left_crop,self.right_crop, self.top_crop,self.bottom_crop
            plt.vlines(x0+lc,y0+bc,y0+tc,color='c',alpha=.7,linewidth=.5)
            plt.vlines(x0+rc,y0+bc,y0+tc,color='c',alpha=.7,linewidth=.5)
            plt.hlines(y0+bc,x0+lc,x0+rc,color='c',alpha=.7,linewidth=.5)
            plt.hlines(y0+tc,x0+lc,x0+rc,color='c',alpha=.7,linewidth=.5,label='Outlier trim')

        ## Draws lines between vials
        if self.debug: print('detector.view_ROI :: Drawing vial/bin lines')
        if bin_lines:
            for item in self.bin_lines[:-1]:
                plt.vlines(self.x + item, y0, y1, color = 'w', alpha = .8, linewidth = 1)
            plt.vlines(self.x + self.bin_lines[-1], y0, y1, color = 'w', alpha = .8, 
                        linewidth = 1, label='Vial boundaries')
        plt.legend()
        plt.tight_layout()
        return

    def display_images(self, cropped_converted, background, subtracted, frame=0,**kwargs):
        '''Generates a three-part figure for manipulated video frames:
          1. Cropped, converted, and grayscaled frame
          2. Null background (median pixel intensity across all indicated (blank_) frames)
          3. Result of subplots 1 - 2 (background subtracted frame)
        ----
        Inputs:
          cropped_converted (nd-array): Cropped and converted nd-array
          background (array): Null background array
          subtracted (nd-array): Background subtracted nd-array
          frame (int): Specific frame/slice of cropped_converted and subtracted
          **kwargs: Arguments for plt.imshow
        ----
        Returns:
          None -- Generates a plot saved in step_3()
          '''
        if self.debug: print('detector.displaying_images: ',end='')
        plt.figure(figsize = (6,8))
    
        ## Displaying the test frame image
        plt.subplot(311)
        if self.debug: print('| Cropped and converted |',end='')
        plt.title('Cropped and converted, frame: %s' % str(frame))
        plt.imshow(cropped_converted[frame], cmap = cm.Greys_r, **kwargs)
        plt.ylabel('Pixels')

        ## Displaying the background image
        plt.subplot(312)
        if self.debug: print(' Background image |',end='')
        plt.title('Background image')
        plt.imshow(background, cmap = cm.Greys_r, **kwargs)
        plt.ylabel('Pixels')

        ## Displaying the background subtracted, test frame image
        plt.subplot(313)
        if self.debug: print(' Subtracted background')
        plt.title('Subtracted background')
        plt.imshow(subtracted[frame], cmap = cm.Greys_r, **kwargs)
        plt.xlabel('Pixels')
        plt.ylabel('Pixels')

        plt.tight_layout()
        return

    def image_metrics(self, spots, image, metric, colorbar=False, **kwargs):
        '''Creates a plot with spot metrics placed over the video image
        ----
        Inputs:
          spots (DataFrame): DataFrame (df_big) containing the spots and their metrics
          image (array): Image background for scatter point plot
          metric (str): Spot metric to filter for
          colorbar (bool): Include a color bar legend
          **kwargs: Keyword arguments to use with plt.scatter
        ----
        Returns:
          None'''
        if self.debug: print('detector.image_metrics')

        ## Create plot
        plt.title(metric)
        plt.imshow(image, cmap = cm.Greys_r)
        plt.scatter(spots.x,spots.y, c = spots[metric], cmap = cm.coolwarm, **kwargs)
        
        ## Add in colorbar
        if colorbar: plt.colorbar()
        
        ## Format plot
        plt.ylim(self.h,0)
        plt.xlim(0,self.w)
        plt.tight_layout()
        return

    def colored_hist(self, spots, metric, bins=40, predict_threshold=False, threshold=None):
        '''Creates a colored histogram to go with image_metrics.
        ----
        Inputs:
          spots (DataFrame): DataFrame with spot metrics (df_big)
          metric (str): Spot metric to evaluate
          bins (int): Number of bins to include for histogram. Default = 40
          predict_threshold (bool): Will predict a threshold for 'signal' metric
          threshold (int): Filtering threshold
        ----
        Returns:
          None'''
        if self.debug: print('detector.colored_hist')        
        
        ## Testing threshold input value
        try: threshold = int(threshold)
        except: pass
    
        ## Assembling histogram parameters
        set_cm = plt.cm.get_cmap('coolwarm')
        n, bin_assignments, patches = plt.hist(spots[metric],bins = bins)
        bin_centers = 0.5 * (bin_assignments[:-1] + bin_assignments[1:])
        col = bin_centers - min(bin_centers)
        col /= max(col)

        ## Plotting by color
        for c, p in zip(col, patches):
            plt.setp(p, 'facecolor', set_cm(c))
    
        ## Getting height of vertical line
        y_max = np.histogram(spots[metric],bins=bins)[0].max()    

        ## Plotting vertical line for eccentricity and mass
        if metric == 'ecc': x_pos = [self.ecc_low,self.ecc_high]
        if metric == 'mass': x_pos = [self.minmass]
        if metric == 'ecc' or metric == 'mass': plt.vlines(x_pos,0,y_max,color = 'gray')
            
        ## Estimate auto-threshold
        if predict_threshold:
            _threshold = self.find_threshold(spots[metric],bins = bins)
            plt.vlines(_threshold,0,y_max,color = 'gray',label='Auto')

        ## Add in user-defined threshold vs. auto
        if isinstance(threshold, int) or isinstance(threshold, float):
            plt.vlines(threshold,0,y_max,color = 'k', label='User-defined')
            if predict_threshold:
                plt.legend(frameon=False)
        
        ## Adding y-axis labels
        plt.ylabel("Counts")
        return

    def spot_checker(self, spots, metrics=['signal'], image=None, **kwargs):
        '''Generates figure containing subplots for image_metrics and colored_hist for
          different spot metrics
        ----
        Inputs:
          spots (DataFrame): DataFrame with spot metrics (df_big)
          metrics (list): List of the metrics to include when generating the plots
          image(array): Background image for plots, default = clean_stack[0]
          **kwargs: Keyword arguments to use with plt.imshow in image_metrics function
        ----
        Returns:
          None'''
        if self.debug: print('detector.spot_checker')

        ## Setting up figure parameters
        subplot = int(str(len(metrics)) + str(2) + str(0))
        if image==None: image = self.clean_stack[0]
        count = 0
        plt.figure(figsize=(4+image.shape[1]/150,len(metrics)*2))
        
        ## Plotting each of the detector spot results over the image
        for i in range(0,len(metrics)):
            ## Defining spot metric and whether to auto-threshold
            col = metrics[i]
            if col=='signal': predict = True
            else: predict = False

            ## Drawing histogram, showing distribution of spots by metric
            count += 1
            plt.subplot(len(metrics),2,count)
            plt.title('Histogram for: %s' % col)
            plt.xlabel(col)
            self.colored_hist(spots,metric=col, bins = 40,predict_threshold=predict)
    
            ## Drawing image plot, colored by metric
            count += 1
            plt.subplot(len(metrics),2,count)
            plt.title('Spot overlay: %s' % col)
            self.image_metrics(spots,image, metric=col,**kwargs)
            plt.ylabel('Pixels')
            if col=='signal':
                plt.xlabel('Pixels')

        plt.tight_layout()
        return

    def find_spots(self, stack, diameter = 3, quiet=True,**kwargs):
        '''Locates the x,y-coordinates for all spots across frames
        ----
        Inputs:
          stack (nd-array): cropped, grayscaled, and background subtracted nd-array
          diameter (int): Estimated diameter of a spot, odds only
          quiet (bool): True silences the output
          **kwargs: Keyword arguments to use with trackpy.batch
        ----
        Returns:
          spots (DataFrame): DataFrame containing all the spots from the TrackPy output
        '''
        if self.debug: print('detector.find_spots')
        ## Check diameter
        diameter = int(diameter)
        if ~diameter % 2: diameter = diameter + 1
    
        ## Option to silence output
        if quiet: tp.quiet()
    
        ## Detect spots
        spots = tp.batch(stack,diameter = diameter, **kwargs)
        
        ## Sorting DataFrame
        spots = spots[spots.raw_mass > 0].sort_values(by='frame')
        return spots
                
    def particle_finder(self, invert=True, **kwargs):
        '''Finds spots and formats the resulting DataFrame. Output can be used with TrackPy.
        ----
        Inputs:
          invert (bool): True if light background, False if dark background
          **kwargs: Keyword arguments to use with trackpy.batch
        ----
        Returns:
          df (DataFrame): DataFrame containing spots and their metrics, becomes df_big'''
        if self.debug: print('detector.particle_finder')

        ## Main spot detection function
        df = self.find_spots(stack = self.spot_stack,
                             quiet=True,invert=True,
                             diameter=self.diameter,
                             minmass=self.minmass,
                             maxsize=self.maxsize)

        ## Catch if there are no spots detected
        if df.shape[0] == 0:
            print('!! Skipping video: No spots detected. Try modifying diameter, maxsize, or minmass parameters')
            raise SystemExit
        
        ## Rounding detector outputs
        df['x'] = df.x.round(2)
        df['y'] = df.y.round(2)
        df['t'] = [round(item/self.frame_rate,3) for item in df.frame]
        df['mass'] = df['mass'].astype(int)
        df['size'] = df.size.round(3)
        df['ecc'] = df.ecc.round(3)
        df['signal'] = df.signal.round(2)
        df['raw_mass'] = df['mass'].astype(int)
        df['ep'] = df.ep.round(1)
        df['True_particle'] = np.repeat(True,df.shape[0])
        return df

    def find_threshold(self,x_array,bins=40):
        '''Auto-generates a signal threshold by finding the local minimum between two local
          maxima, or takes the average between 0 and the global maximum
        ----
        Inputs:
          x_array (list): list of all metric (signal) points to find a threshold for
          bins (int): Number of bins to search across, default = 40.
        ----
        Returns:
          threshold (int): Auto-generated threshold
        '''
        if self.debug: print('detector.find_threshold')
        
        ## Looking at the distribution of spot metrics as a histogram
        x_array = np.histogram(x_array,bins = bins)[0]
        
        ## Peak finding with SciPy.signal module
        peaks = find_peaks(x_array)[0]
        prominences = peak_prominences(x_array,peaks)
        threshold = find_peaks(x_array, prominence=np.max(prominences))[0][0]
        if self.debug: print('                   Threshold =',threshold)
        return threshold

    def invert_y(self,spots):
        '''Inverts spots along the y-axis. Important for converting spots indexed for an image to a plot.
        ----
        Inputs:
          spots (DataFrame): DataFrame containing a 'y' column
        ----
        Returns:
          inv_y (list): In-place list of the y-coordinates inverted'''
        if self.debug: print('detector.invert_y')
        
        ## Inverts y-axis
        inv_y = abs(spots.y - spots.y.max())
        return inv_y
    
    def get_slopes(self):
        '''Creates a dictionary with keys for vials and values for the DataFrame sliced by
        vial. It will also calculate the local linear regression for each vial and
        returns the DataFrame containing all the slopes and linear regression statistics
        for that vial.
        ----
        Inputs (Imported from the detector object):
          df (DataFrame) : DataFrame sliced by vial
          vials (int) : Number of vials in video
          window (int) : Window size to calculate local linear regression
          vial_ID (str) : Vial-specific ID, taken from the first 'vial_id_vars' of naming convention
        ----
        Returns (Exported tp the detector object):
          result (dict) : DataFrame containing the local linear regression statistics
            by vial
          vial (dict) : Dictionary of DataFrames sliced by vial, keys are vials and values
            are DataFrames'''
        if self.debug: print('detector.get_slopes')

        ## Create empty dictionaries
        self.vial,self.result = dict(),dict()

        ## Slicing DataFrame (df.filtered) by vial and assigning slices to dictionary keys (vials)
        for i in range(1,self.vials + 2):
            try:
                ## Set dict key to '1' if only 1 vial, otherwise set dict key to vial number
                if self.vials == 1 or i == self.vials + 1: self.vial[i] = self.df_filtered
                else: self.vial[i] = self.df_filtered[self.df_filtered.vial==i]
            
                ## Setting the result to the result from the local linear regression
                self.result[i] = self.local_linear_regression(self.vial[i]).iloc[0].tolist()
            
                ## Add vial_ID to the result
                if i == self.vials + 1: v = 'all'
                else: v = i
                
                ## Name vial_ID
                vial_ID = ['_'.join(self.vial_ID) + '_'+ str(v)]
                self.result[i] = vial_ID + self.result[i]
                
                ## Rounding results so they are more manageable and require less space.
                self.result[i][1:3] = [int(item) for item in self.result[i][1:3]]
                self.result[i][3:] = [round(item,4) for item in self.result[i][3:]]
            except:
                print('Warning:: Could not process vial %s' % i)
        return

    def get_trim_lines(self,df,edge = 'top',sensitivity=1):
        '''Calculates spacial thresholds for cropping outlier points at the edge of window
        ----
        Inputs:
          df (DataFrame): DataFrame of all points to consider
          edge (str) {'top'|'bottom','left','right'}: Which edge to trim
           sensitivity (float): How sensitive to make the thresholding
        ----
        Returns:
          crop (float): Cutoff threshold'''
        if self.debug: print('detector.get_trim_lines ::')
        for _edge in edge:
            _list,diff_list = [],[]

            # Define axis
            if edge == 'top' or edge == 'bottom': axis = 'y'
            if edge == 'left' or edge == 'right': axis = 'x'

            # Define quantile starting value
            if edge == 'left' or edge == 'bottom': quant = 0
            if edge == 'right' or edge == 'top': quant = 0.96

            ## Find quantile boundaries
            for i in range(5):
                val = df[axis].quantile(quant + i * 0.01)
                _list.append(val)

            ## Get difference between quantile boundaries
            for i in range(4):
                diff_list.append(abs(_list[i]-_list[i+1]))

            ## Calculate boundary, cutoff, and thresholds
            ## --boundary as most extreme value
            ## --cutoff as median 0-4 or 95-99 percentiles x scalar (sensitivity)
            ## --threshold as difference between boundary and cutoff
            ## --crop as value to crop points at 
        
            if 'right' in edge or 'top' in edge:
                boundary = _list[-1] # Max value
                cutoff = np.median(diff_list[:-1]) * sensitivity
                threshold = boundary - cutoff

                if cutoff < threshold: crop = boundary - cutoff
                else: crop = cutoff
                if self.debug: print('detector.get_trim_lines ::',edge,'@',crop)
                return crop

            if 'left' in edge or 'bottom' in edge:
                boundary = _list[0] # Min value
                cutoff = np.median(diff_list[1:]) * sensitivity 
                threshold = boundary + cutoff
            
                if cutoff > threshold: crop = boundary + cutoff
                else: crop = boundary
                if self.debug: print('detector.get_trim_lines ::',edge,'@',crop,'(no crop)')            
                return crop
                
    def bin_vials(self, df, vials, percentage=1,top=False, bin_lines=None):
        '''Bin spots into vials. Function takes into account all points along the x-axis, 
          and divides them into specified number of bins based on the min and max
          points in the array.
        ----
        Inputs:
          df (DataFrame):
          vials (int): Number of vials in video
        Returns:
          bin_lines (list): Binning intervals along the x-axis
          spot_assignments (pd.Series): '''
        if self.debug: print('detector.bin_vials')

        ## Bin vials, conditional for vial quantity
        if vials == 1:
            if type(bin_lines) == 'list': bin_lines = bin_lines
            else: bin_lines = [df.x.min(),df.x.max()] 
            spot_assignments = np.repeat(1,df.shape[0])
        else: ## More than 1 vial
            if type(bin_lines) == 'list': bin_lines = bin_lines
            else: bin_lines = pd.cut(df.x,vials,include_lowest=True,retbins=True)[1]

            ## Assign spots to vials
            _labels = range(1,vials+1)
            spot_assignments = pd.cut(df.x, bins=bin_lines, labels=_labels)
            spot_assignments = pd.Series(spot_assignments).astype('int')
        
            ## Checks to make sure all vials have at least one spot. Important if a middle vial is absent or vials binned incorrectly.
            counts = np.unique(spot_assignments, return_counts = True)
            for v,c in zip(counts[0], counts[1]):
                if c == 0:
                    print('Warning: vial',v,'is empty and cannot be evaluated')

        return bin_lines, spot_assignments


    def local_linear_regression(self, df, method = 'max_r'):
        '''Performs a local linear regression, using a user-defined sliding window (self.window)
        ----
        Inputs:
          df (DataFrame): DataFrame containing all formatted points (df_filtered). 'y' needs to be converted from image indexing to plot indexing
          method (str): Two options: greatest regression coefficient (max_r) or lowest error (min_err)
        ----
        Returns:
          result (DataFrame): Single-row slice of a DataFrame corresponding with the results from the local linear regression
        '''
        if self.debug: print('detector.local_linear_regression')

        ## Defining empty variables
        result_list, result = [],pd.DataFrame()
        llr_columns = ['first_frame','last_frame','slope','intercept','r','pval','err']#,'count_llr','count_all']
        
        _count_all = np.median(df.groupby('frame').frame.count())
        
        ## Iterating through the window
        frames = (self.crop_n - self.crop_0) - self.window
        for i in range(frames):
            ## Defining search parameters for each iteration
            start, stop = int(i),int(i+self.window)
            df_window  = df[(df.frame >= start) & (df.frame <= stop)]

            ## Testing if there are enough frames in the slice
            if df_window.groupby('frame').y.count().min() == 0:
                print('Issue with number of frames with flies detected vs. window size')
                print(i, i + self.window, len(df_window.frame.unique()))
                continue
                        
            ## Performing local linear regression
            try: 
                ## Grouping points by frame
                _frame = df_window.groupby('frame').frame.mean()
                _pos  = df_window.groupby('frame').y.mean()
                _count_llr = np.median(df_window.groupby('frame').frame.count())

                ## Performing linear regression on subset and formatting output to list
                _result = linregress(_frame,_pos)
                _result = [start,stop] + np.hstack(_result).tolist() #+ [_count_llr,_count_all]

                ## If slope is not significantly different from 0, then set slope = 0
                if _result[-2] >= 0.05: _result[2] = 0

            ## Have row of NaN if unable to process
            except: _result = [start,stop] + [np.nan,np.nan,np.nan,np.nan]
#             except: _result = [start,stop] + [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]

            ## Add results list to a list of lists
            result_list.append(_result)
        
        ## Assembles the list of lists into a DataFrame
        result = pd.DataFrame(data=result_list,columns=llr_columns)

        ## Filtering method
        if method == 'max_r': result = result[result.r == result.r.max()]
        elif method == 'min_err': result = result[result.err == result.err.min()]
        else: print("Unrecognized method, chose either 'max_r' to select the window with the greatest R or 'min_err' to select the window with the lowest error")
        return result
        
        
    def step_1(self, gui = False, grayscale = True):
        '''Crops and formats the video, previously loaded during detector initialization.
        ----
        Inputs:
          gui (bool): True creates plots for detector optimization
          grayscale (bool): True converts video array to grayscale
        ----
        Returns:
          None'''
        print('-- [ Step 1  ] Cleaning and format image stack')
        x,y = self.x,self.y
        x_max, y_max = int(x + self.w),int(y + self.h)
        stack = self.image_stack
        
        if self.debug: print('detector.step_1 cropped and grayscale: grayscale image:', grayscale)
        self.clean_stack = self.crop_and_grayscale(stack,
                     y=y, y_max=y_max,
                     x=x, x_max=x_max,
                     first_frame=self.crop_0, 
                     last_frame=self.crop_n,
                     grayscale=grayscale)

        ## Confirm frame ranges
        self.check_variable_formats()
        if self.blank_0 < self.crop_0:
            self.blank_0 = self.crop_0
        if self.blank_n > self.crop_n:
            self.blank_n = self.crop_n
        
        if grayscale:
            if self.debug: print('detector.step_1 cropped and grayscale: grayscale image')
            self.clean_stack = self.crop_and_grayscale(stack,
                         y=y, y_max=y_max,
                         x=x, x_max=x_max,
                         first_frame=self.crop_0, last_frame=self.crop_n)
        else:
            if self.debug: print('detector.step_1 cropped and grayscale: no color image')
            self.clean_stack = self.crop_and_grayscale(stack,
                         y=y, y_max=y_max,
                         x=x, x_max=x_max,
                         first_frame=self.crop_0, last_frame=self.crop_n, grayscale=False)                        

        if self.debug: print('detector.step_1 cropped and grayscale dimensions: ', self.clean_stack.shape)

        ## Subtracts background to generate null background image and spot stack
        self.spot_stack,self.background = self.subtract_background(video_array=self.clean_stack)
        if self.debug: print('detector.step_1 spot_stack and null background created')
        return


    def step_2(self):
        '''Performs spot detection and manipulates the resulting DataFrames'''
        print('-- [ Step 2  ] Identifying spots')

        ## Particle detection step
        self.df_big = self.particle_finder(minmass=self.minmass,diameter=self.diameter,
                                            maxsize=self.maxsize, invert=True)
        if self.debug: print('                   Identified %s spots' % self.df_big.shape[0])
        return


    def step_3(self, gui = False):
        '''Visualizes spot metrics
        ----
        Inputs:
          gui (bool): True creates plots for detector optimization
        ----
        Returns:
          None
        '''
        print('-- [ Step 3  ] Visualize spot metrics ::',gui)
        if gui:        
            ## Visualizes spot metrics on plot with accompanying color-coded histogram
            self.spot_checker(self.df_big,metrics=['ecc','mass','signal'], alpha=.1)
            
            plot_spot_check = self.name_nosuffix + '.spot_check.png'
            plt.savefig(plot_spot_check,dpi=200)
            print('                --> Saved:',plot_spot_check.split('/')[-1])
            plt.close()
        
            ## Creating image plot with rectangle superimposed over first frame
            plt.figure()
            self.display_images(self.clean_stack,self.background,self.spot_stack,frame=20)
            plt.tight_layout()
            plot_name = self.name_nosuffix + '.processed.png'
            plt.savefig(plot_name, dpi=100)
            plt.close()
            print('                --> Saved:',plot_name.split('/')[-1])
        return

    def step_4(self):
        '''Filters and processes data detected points'''

        ## Assigning spots a True/False status based on ecc/eccentricity (circularity)
        def ecc_filter(x,low=0,high=1):
            '''Simple function for making a vector True/False depending on an upper and lower bound
            ----
            Inputs:
              x (numeric): value to perform operation on
              low (numeric): Lower bound
              high (numeric): Upper bound
            Returns:
              (bool): True/False'''
            if x >= low and x <= high: return True
            else: return False

        print('-- [ Step 4a ]   - Setting spot threshold')        
        ## Auto-detecting threshold
        if self.threshold == 'auto': self.threshold = self.find_threshold(self.df_big.signal)

        print('-- [ Step 4b ]   - Filtering by signal threshold') 
        ## Assigning spots a True/False status based on signal threshold
        self.df_big['True_particle'] = [x >= self.threshold for x in self.df_big.signal]

        t_or_f = np.unique(self.df_big.True_particle, return_counts=True)
        if self.debug: print('                   True (%s) and False (%s) spots' % (t_or_f[0],t_or_f[1]))
        
        print('-- [ Step 4c ]   - Filtering by eccentricity/circularity') 
        self.df_big.loc[self.df_big.True_particle == True,'True_particle'] = self.df_big[self.df_big.True_particle==True].ecc.map(lambda x: ecc_filter(x,low=self.ecc_low,high=self.ecc_high))
        t_or_f = np.unique(self.df_big.True_particle, return_counts=True)
        if self.debug: print('                   True (%s) and False (%s) spots'%(t_or_f[0],t_or_f[1]))
        
        ## Checking to confirm DataFrame is not empty after filtering
        if self.df_big[self.df_big.True_particle].shape[0] == 0:
            print('\n\n!! No spots post-filtering, check detector and background subtraction settings for proper optimization')
            raise SystemExit

        ## Pruning errant points on periphery if outliers
        print('-- [ Step 4d ]   - Trimming outliers (if indicated)')
        if self.trim_outliers:
            self.left_crop = self.get_trim_lines(self.df_big,edge='left',sensitivity = self.outlier_LR)
            self.right_crop = self.get_trim_lines(self.df_big,edge='right',sensitivity = self.outlier_LR)
            self.top_crop = self.get_trim_lines(self.df_big,edge='top',sensitivity = self.outlier_TB)
            self.bottom_crop = self.get_trim_lines(self.df_big,edge='bottom',sensitivity = self.outlier_TB)
            
            self.df_big = self.df_big[(self.df_big.x >= self.left_crop) & (self.df_big.x <= self.right_crop) &
                                        (self.df_big.y <= self.top_crop) & (self.df_big.y >= self.bottom_crop)]
        
        ## Assigning spots to vials, 0 if False AND outside of the True point range
        print('-- [ Step 4e ]   - Assigning spots to vials')
        self.bin_lines, self.df_big.loc[self.df_big['True_particle'],'vial'] = self.bin_vials(self.df_big[self.df_big.True_particle],vials = self.vials)
        
        ########################################
        ## Beginning of publication insert
        if publication:
            df=self.df_big
            self.bin_lines = self.bin_vials(df[df.y > 120], vials= self.vials)[0]
            df['vial'] = np.repeat(0,df.shape[0])
            vial_assignments = self.bin_vials(df, vials = self.vials, bin_lines = self.bin_lines)[1]
            df.loc[(df.x >= self.bin_lines[0]) & (df.x <= self.bin_lines[-1]),'vial'] = vial_assignments
            self.df_big = df
        ## End of publication insert
        ########################################
        
        self.df_big.loc[self.df_big['True_particle']==False,'vial'] = 0

        print('-- [ Step 4f ]   - Saving raw data file')
        ## Saving the TrackPy results, plus filter and vial notations        
        self.df_big.to_csv(self.path_data, index=None)
        print('                --> Saved:',self.path_data.split('/')[-1])

        return

    def step_5(self):
        '''Calculates local linear regressions'''
        print('-- [ step 5  ] Setting up DataFrames for local linear regression')

        ## Filtering spots and pruning unnecessary columns
        self.df_filtered = self.df_big[(self.df_big.True_particle) & (self.df_big.vial != 0)]
        if self.debug: print('self.df_filtered.shape:',self.df_filtered.shape)
        self.df_filtered = self.df_filtered.drop(['ecc','signal','ep','raw_mass','mass','size','True_particle'],axis=1)
        self.df_filtered = self.df_filtered.sort_values(by=['vial','frame','y','x'])    

        ## Adding experimental details to DataFrame
        self.specify_paths_details(self.video_file)
       
       ## Filling in experimental details to DataFrame
        for item in self.file_details.keys():
            self.df_filtered[item] = np.repeat(self.file_details[item],self.df_filtered.shape[0])        
        
        ## Invert y-axis -- images indexed upper left to lower right but converting because plots got left left to upper right
        self.df_filtered['y'] = self.invert_y(self.df_filtered)
        self.df_filtered['y'] = self.df_filtered.y.round(2)
        
        ## Convert vial assignments from float to int
        self.df_filtered['vial'] = self.df_filtered['vial'].astype('int')

        ## Save the filtered DataFrame
        path_filtered = self.name_nosuffix+'.filtered.csv'
        self.df_filtered.to_csv(self.path_filtered, index=False)
        print('                --> Saved:',self.path_filtered.split('/')[-1])
        return
        
    def step_6(self,gui=False):
        '''Creating diagnostic plots to visualize spots at beginning & end of most linear 
          section, throughout the video, and a vertical velocity plot for each vial.
        ----
        Inputs:
          gui (bool): True creates plots for detector optimization
        ----
        Returns:
          None'''
        print('-- [ step 6a ] Visualize spot metrics ::',gui)
        
        ## Check window size is not greater than video length
        video_length = self.crop_n - self.crop_0
        if self.window > video_length:
            print('!! Issue with window size > video length: was %s, now %s' % (self.window, video_length-1))
            self.window = video_length - 1
            
        if gui:        

            ## Visualizes the region of interest
            plt.figure()
            self.view_ROI(border = True,
                            x0 = self.x, x1 = self.x + self.w,
                            y0 = self.y, y1 = self.y + self.h,
                            bin_lines = True)
        
            plot_roi = self.name_nosuffix + '.ROI.png'
            plt.savefig(plot_roi,dpi=100)
            print('                --> Saved:',plot_roi.split('/')[-1])
            plt.close()
            
        print('-- [ step 6b ] Creating diagnostic plot file')
        ## Set up plots
        plt.figure(figsize=(10,8))
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

        ## Finding the frames that flank the most linear portion 
        ##    of the y vs. t curve for all points, not just by vials
        if self.debug: print('-- [ step 6b ] Plotting data: Re-running local linear regression on all')
        _result = self.local_linear_regression(self.df_filtered)
        begin = _result.iloc[0].first_frame.astype(int)
        end = _result.iloc[0].last_frame.astype(int)
        
        ## For future release
#         min_R = _result.iloc[0].r_value ##

        ## Need only True spots, but not inverted -- df_big
        spots = self.df_big[self.df_big.True_particle]

        ## Creating the diagnostic plot
        print('-- [ step 6b1] Plotting image plots with overlaying points')        
        if self.debug: print("-- [ step 6b1] Plotting data: Plot 1 - Frame %s" % begin)
        self.image_plot(df = spots,frame = begin, ax=ax1)

        if self.debug: print("-- [ step 6b1] Plotting data: Plot 2 - Frame %s" % end)
        self.image_plot(df = spots, ax=ax2, frame = int(str(end)))

        if self.debug: print("-- [ step 6b1] Plotting data: Plot 3 - Frame ALL")
        self.image_plot(df = spots,ax=ax4, frame = None)

        print('-- [ step 6b2] Performing local linear regression')
        self.get_slopes()
        
        print('-- [ step 6b3] Plotting local linear regression results')
        self.loclin_plot(ax=ax3)

        ## Saving diagnostic plot
        fig.tight_layout()
        plt.savefig(self.path_diagnostic,dpi=300, transparent = True)
#         plt.savefig(self.path_project + '/diagnostic_plots/' + self.name + '.diagnostic.png',dpi=100, transparent = True)        

        ## Future release
#         if min_R < self.review_R:
#             plt.savefig(self.path_review_diagnostic,dpi=100, transparent = True)
            
        plt.close()
        plt.close()
        print('                --> Saved:',self.path_diagnostic.split('/')[-1])
        return
        
    def step_7(self):
        '''Writing the video's slope file'''
        print('-- [ step 7  ] Setting up slopes file')
        slope_columns = ['vial_ID','first_frame','last_frame','slope','intercept','r_value','p_value','std_err']#,'count_llr','count_all']
        
        ## Converting dictionary of local linear regressions into a DataFrame
        self.df_slopes = pd.DataFrame.from_dict(self.result,orient='index',columns = slope_columns)

        ## Applying conversion factor if indicated, if not it will just be '1'
        self.df_slopes['slope'] = self.df_slopes.slope.transform(lambda x: x * (self.conversion_factor)).round(4)
        
        ## Adding in experimental details from naming convention into the slopes DataFrame
        for item in self.file_details.keys():
            self.df_slopes[item] = np.repeat(self.file_details[item],self.df_slopes.shape[0])

        ## Specifying column names
        slope_columns = [item for item in self.file_details.keys()] + slope_columns
        self.df_slopes = self.df_slopes[slope_columns]
    
        ## Saving slope file
        self.df_slopes.to_csv(self.path_slope,index=False)
        print('                --> Saved: %s \n' % self.path_slope.split('/')[-1])
        plt.close('all')
        
        print(self.df_slopes[['vial_ID','slope','r_value']])
        print('\n')
        return
        
    def image_plot(self,df,frame=None,ax=None,ylim=[0,1000]):
        '''Image subplot for the diagnostic plot
        ----
        Inputs:
          df (DataFrame): DataFrame containing all spots
          frame (int): Frame to slice df
          ax (int): plot coordinate
          ylim (2-item list): y-limits
        ----
        Returns:
          ax (object): matplotlib object containing plot'''
        if self.debug: print('detector.image_plot')
        
        ## Get frame number
        try: frame = int(frame)
        except: frame = None

        ## Assign plotting parameters depending on which frame(s)
        if type(frame) == int and frame in df.frame.unique():
            df = df[(df.frame == frame)]
            alpha = .25
            title  = "Frame: %s" % frame
            ax.set_title(title)
        elif frame == None:
            frame = 0
            alpha = 0.01
            title = 'All x,y-points throughout video'
            ax.set_title(title)
        elif type(frame) == int and not frame in df.frame.unique():
            print('Error: input frame is not accounted for in current DataFrame')
        else:
            print("Chose a frame value in integer form, or 'None'")

        ## Issue with plotting if last frame in stack
        if frame == self.n_frames: image = self.clean_stack[frame-1]
        else: image = self.clean_stack[frame]

        ## Plotting image
        ax.imshow(image,cmap=cm.Greys_r,origin='upper')
        ax.set_ylim(self.h,0)
        ax.set_xlim(0,self.w)
        
        ## Plotting vertical bin lines
        ax.vlines(self.bin_lines,0,image.shape[0],alpha = .3)  
        
        ## Coloring spots by vial
        df = df.sort_values(by='vial')
        df = df[df.vial != 0]
        if self.vials >= 1:
            ax.scatter(df.x, df.y, 
                        s = 30, 
                        alpha = alpha,
                        c = df['vial'], 
                        cmap=self.vial_color_map)

        ## Getting rid of axis labels
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        return ax
    
    def loclin_plot(self,ax = None):
        '''Local linear regression plot: mean y-position vs. frame or time. Colored by vial
          and bolded for the most linear section
        ----
        Inputs:
          ax (int): plot coordinate
        ----
        Returns:
          None'''
        if self.debug: print('detector.loclin_plot')
        
        def two_plot(df,vial,label,first,last,ax=None):
            '''Adds the bolded flair to the local linear regression plot'''
            if self.debug: print('detector.two_plot')
            
            ## All points
            x = df.groupby('frame').frame.mean()
            y = df.groupby('frame').y.mean()
            
            ## Convert to cm / sec
            if self.convert_to_cm_sec:
                x = x / self.frame_rate
                y = y / self.pixel_to_cm
            ax.plot(x,y, 
                    alpha = .35, 
                    color = self.color_list[vial-1],
                    label='') 

            ## Only points in the most linear segment
            df = df[(df.frame >= first) & (df.frame <= last)]
            x = df.groupby('frame').frame.mean()
            y = df.groupby('frame').y.mean()
            
            ## Convert to cm / sec
            if self.convert_to_cm_sec:
                x = x / self.frame_rate
                y = y / self.pixel_to_cm
            ax.plot(x,y,
                    color = self.color_list[vial-1],
                    label = label)
            return
    
        ## Plotting multiple vials' data
        if len(self.result) > 1:
            for V in range(1,self.vials+1):
                l = 'Vial '+str(V)
                c = self.color_list[V-1]
                _details = self.result[V]
                frames = _details[1],_details[2]
                two_plot(self.vial[V],
                         vial = V,
                         label = l,
                         first = _details[1],
                         last  = _details[2],
                         ax=ax)

        ## Add on labels and legends
        ax.legend(loc=2, frameon=False, fontsize='x-small')
        label_y,label_x = 'pixels','Frames'
        if self.convert_to_cm_sec: 
            label_x = 'Seconds'
            label_y = 'cm'
        title = 'Cohort climbing kinematics'
        label_y = 'Mean y-position (%s)' % label_y
        ax.set_ylim(0,ax.get_ylim()[1])
        ax.set(title=title,xlabel = label_x,ylabel=label_y)
        
        if ax == None: return
        else: return ax

    ## Parameter testing is only used in the GUI
    def parameter_testing(self, variables, axes):
        '''Parameter testing in the GUI and done separately to account for plots with wx'''
        if self.debug: print('detector.parameter_testing')
        
        ## Running through the first few steps
        self.load_for_gui(variables)
        
        ## Load in video
        self.step_1(gui=True) # Crop and convert video

        ## Detect spots
        self.step_2()
        
        ## First optimization plot
        self.step_3(gui=True)
        
        ## Filters DataFrame of detected spots
        self.step_4() 

        ## Executing final steps
        self.step_5()
        self.step_6(gui=True)
        self.step_7()
        
        #### Working through the GUI plots        
        ## Setting plot (upper left) for background image
        if self.debug: print('detector.parameter_testing: Subplot 0: Background image')
        axes[0].set_title("Background Image")
        axes[0].imshow(self.background,cmap=cm.Greys_r)
        axes[0].set_xlim(0,self.w)
        axes[0].set_ylim(self.h,0)
        axes[0].scatter([0,self.w],[0,self.h],alpha=0,marker='.')
        
        ## Slice df_big into true vs. false spots
        if self.debug: print('detector.parameter_testing: Slicing DataFrames')
        spots_false = self.df_big[~self.df_big['True_particle']]
        spots_true = self.df_big[self.df_big['True_particle']]
        
        ## Binning and coloring spots
#         bin_lines,spots_true['vial'] = self.bin_vials(spots_true,vials = self.vials)
#         spots_true = spots_true[(spots_true.x >= self.bin_lines.min()) & (spots_true.x <= self.bin_lines.max())]

        spots_true['vial'] = np.repeat(0,spots_true.shape[0])
        vial_assignments = self.bin_vials(spots_true, vials = self.vials, bin_lines = self.bin_lines)[1]
        spots_true.loc[(spots_true.x >= self.bin_lines[0]) & (spots_true.x <= self.bin_lines[-1]),'vial'] = vial_assignments

        spots_true.loc[:,'color'] = spots_true.vial.map(dict(zip(range(1,self.vials+1), self.color_list)))
        bins=40

        ## Setting plots for scatterplot overlay on a selected frame
        if self.debug: print('detector.parameter_testing: Subplot 1: Test frame')
        axes[1].set_title('Frame: '+str(self.check_frame))
        axes[1].imshow(self.clean_stack[self.check_frame], cmap = cm.Greys_r)
        axes[1].scatter(spots_false[(spots_false.frame==self.check_frame)].x,
                        spots_false[(spots_false.frame==self.check_frame)].y, 
                        color = 'b',marker ='+',alpha = .5)
        a = axes[1].scatter(spots_true[spots_true.frame==self.check_frame].x,
                            spots_true[spots_true.frame==self.check_frame].y, 
                            c = spots_true[spots_true.frame==self.check_frame].vial,
                            cmap = self.vial_color_map,
                            marker ='o',alpha = .8)
        a.set_facecolor('none')
        axes[1].vlines(self.bin_lines,0,self.df_big.y.max(),color='w')
        axes[1].set_xlim(0,self.w)
        axes[1].set_ylim(self.h,0)
        
        ##########
        ## Fly counts
#         axes[5].plot(spots_true.groupby('frame').frame.mean(),
#                      spots_true.groupby('frame').frame.count(), 
#                      label = 'Fly count', color = 'g',alpha = .5)
#                      
#         axes[5].hlines(np.median(spots_true.groupby('frame').frame.count()),
#                        self.df_big.frame.min(),self.df_big.frame.max(),
#                        linestyle = '--',alpha = .5, 
#                        color = 'gray', label = 'Median no. flies')
#         axes[5].set(title = 'Flies per frame',
#                         ylabel='Flies detected',
#                         xlabel='Frame') 
#         axes[5].legend(frameon=False, fontsize = 'small')

        ##########
        df = self.df_filtered.sort_values(by='frame')
        for V in range(1,self.vials + 1):
            color = self.color_list[V-1]
            _df = df[df.vial == V]

            ## Local linear regression
            begin = self.local_linear_regression(_df).iloc[0].first_frame.astype(int)
            end = begin + self.window        
                        
            ## Plotting all points
            axes[5].plot(_df.groupby('frame').frame.unique(),
                _df.groupby('frame').y.count(),alpha = .3, color = color,label='') 
            
            ## Plotting most linear points
            _df = _df[(_df.frame >= begin) & (_df.frame <= end)]
            axes[5].plot(_df.groupby('frame').frame.unique(),
                _df.groupby('frame').frame.count() ,color = color, alpha = .5)
    
            axes[5].hlines(np.median(_df.groupby('frame').frame.count()),
                       df.frame.min(),df.frame.max(),
                       linestyle = '--',alpha = .7, 
                       color = color)


        # Deciding number of columns for legend
        if self.vials > 10: ncol = 3
        elif self.vials > 5: ncol = 2
        else: ncol=1
        
        ## Setting labels
        label_y,label_x = '(pixels)','Frames'
        if self.convert_to_cm_sec: 
            label_x,label_y = 'Seconds','(cm)'
        labels = ['Flies detected per frame','Flies detected','Frame']
        axes[5].set(title = labels[0], ylabel=labels[1],xlabel=labels[2]) 
        axes[5].set_ylim(ymin = 0,ymax = np.max(_df.groupby('frame').frame.count())*1.2)
#         axes[5].legend(frameon=False, fontsize='x-small', ncol=ncol)

        custom_lines = [Line2D([0], [0], color='k', linestyle = '--', alpha = .9),
                        Line2D([0], [0], color='k', linestyle = '-', alpha = .5)]
        custom_labels = ['Median', 'All frames']
        axes[5].legend(custom_lines, custom_labels,frameon=False, fontsize='x-small', ncol=ncol)

        #############



        ## Mass histogram
        axes[3].set_title('Mass Distribution')
        axes[3].hist(self.df_big.mass,bins = bins)
        y_max = np.histogram(self.df_big.mass,bins=bins)[0].max()
        axes[3].vlines(self.minmass,0,y_max)

        ## Signal histogram
        axes[4].set_title('Signal Distribution')
        axes[4].hist(self.df_big.signal,bins = bins)
        y_max = np.histogram(self.df_big.signal,bins=bins)[0].max()
        axes[4].vlines(self.threshold,0,y_max)

        ## Calculating local linear regression
#         self.step_5()
        
        ## Setting plots for local linear regression
        df = self.df_filtered.sort_values(by='frame')

        ## Converting to cm per sec if specified
        convert_x,convert_y = 1,1
        if self.convert_to_cm_sec:
            convert_x,convert_y = self.frame_rate,self.pixel_to_cm

        ## LocLin plot for each vial
        for V in range(1,self.vials + 1):
            label = 'Vial '+str(V)
            color = self.color_list[V-1]
            _df = df[df.vial == V]

            ## Local linear regression
            begin = self.local_linear_regression(_df).iloc[0].first_frame.astype(int)
            end = begin + self.window        
            
            ## Plotting all points
            axes[2].plot(_df.groupby('frame').frame.mean() / convert_x,
               _df.groupby('frame').y.mean() / convert_y,alpha = .35, color = color,label='') 
            
            ## Plotting most linear points
            _df = _df[(_df.frame >= begin) & (_df.frame <= end)]
            axes[2].plot(_df.groupby('frame').frame.mean() / convert_x,
           _df.groupby('frame').y.mean() / convert_y,color = color, label = label)

        # Deciding number of columns for legend
        if self.vials > 10: ncol = 3
        elif self.vials > 5: ncol = 2
        else: ncol=1
        
        ## Setting labels
        label_y,label_x = '(pixels)','Frames'
        if self.convert_to_cm_sec: 
            label_x,label_y = 'Seconds','(cm)'
        labels = ['Mean vertical position over time','Mean y-position %s' % label_y,label_x]
        axes[2].set(title = labels[0], ylabel=labels[1],xlabel=labels[2]) 
        axes[2].legend(frameon=False, fontsize='x-small', ncol=ncol)   

        return