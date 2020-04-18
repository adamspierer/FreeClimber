#!/usr/bin/env python
# -*- coding: utf-8 -*- 
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
# from matplotlib.cm import Greys_r, coolwarm,viridis
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap

# import warnings
# warnings.simplefilter(action='ignore', category=FutureWarning)


### To do
## [ ] ArgParse: Option to do all or undone
## [ ] Input/output for lists of files processed/with errors
## [X] Close plot window after saving



class detector(object):
    '''Particle detection platform for identifying the group climbing velocity of a 
    group of flies (or particles) in a Drosophila negative geotaxis (climbing) assay.
    This platform is designed to take videos of varying background homogeneity and 
    extract the x,y,time coordinates of flies by subtracting the static portion of the 
    video. Climbing velocities are calculated as the most linear portion of a user-defined
    subset of frames--done so by vial if defined. 
    '''
    def __init__(self, video_file, config_file = None,gui=False, variables = None):
        '''Initializing detector object'''
        print('---------------> Initializing detector')
        self.video_file = self.check_video(video_file)
        
        if gui:
            self.load_for_gui(variables = variables)
        elif  not gui:
            self.load_for_main(config_file = config_file)

#         for item in variables:
#             print(item)

        print('\n\n')
        self.specify_paths_details(video_file)
        self.image_stack = self.video_to_array(video_file,loglevel='panic')
        return


## Loading functions    
    def load_for_gui(self,variables):
        print('---------------> GUI: loading variables for a GUI run')
        self.config = None
        self.path_project = None

        if variables == None:
            print('\n\nExiting program. No variable list loaded')
            raise SystemExit
        for item in variables:
#             print(item)
            try: exec('self.'+item)
            except: pass
        return     
        
    def load_for_main(self, config_file = None):
        print('---------------> MAIN: loading variables for command line run')
        if config_file == None:
            self.config_file = sys.argv[1]
        self.check_cfg(self.config_file)
        self.load_from_cfg(self.config_file)
        return
        
    def load_from_cfg(self,config_file=None):
        print('---------------> load_from_cfg: Reading in variables from .cfg file')
        try:
            with open(config_file,'r') as f:
                variables = f.read().split('\n')
            f.close()
            
            variables = [item for item in variables if not item.startswith(('#',' ','\n'))]
            for item in variables:
                try: exec('self.'+item)
                except: pass
        except: 
            print('\n\nExiting program. Could not read in file.cfg, but path and suffix are good')
            raise SystemExit
        return

## Checking functions
    def check_video(self,video_file=None):
        print('---------------> Confirming video path is a file')
        if os.path.isfile(video_file): pass
        else:
            print('\n\nExiting program. Invalid path to video file: ',video_file)
            raise SystemExit
        return video_file
        
    def check_cfg(self,config_file=None):
        print('---------------> Check_cfg: confirming .cfg file is a valid path')
        if os.path.isfile(config_file) and config_file.endswith('.cfg'): pass
        else:
            print('\n\nExiting program. Not a valid configuration file')
            raise SystemExit
        return

## Specifying variables
    def specify_paths_details(self,video_file):
        print('---------------> specify_paths_details: creating variables from those already loaded')
        folder,name = os.path.split(video_file)
        self.name = name[:-5]
        self.name_nosuffix = '.'.join(video_file.split('.')[:-1])
        self.experiment_details = self.name.split('_')
        self.file_details = dict(zip(self.naming_convention.split('_'),self.name.split('_')))
        
        ## Defining final file names and destinations
        folder_names = ['data','filter','plot','slope']
        folder_suffixes = ['.raw.csv','.filter.csv','.diag.png','.slopes.csv']
        for item,jtem in zip(folder_names,folder_suffixes):
            var_name = 'self.path_'+item
            file_path = os.path.join(folder,name+jtem)
            exec(var_name+"='"+file_path+"'")
        
        ## Project folder specific paths
        if self.path_project == None: self.path_project = os.path.join(folder,self.name + '.cfg')
#         self.path_result  = os.path.join(self.path_project,'results.csv')

        ## Naming individual vials
        self.experiment_details[-1] = '.'.join(self.experiment_details[-1].split('.')[:-1])
        self.vial_ID = self.experiment_details[:self.vvars]
        self.discrete_cmap(vials=self.vials)
#         print(self.color_list)
        return

## Video processing functions
    def video_to_array(self, file, **kwargs):
        print('---------------> video_to_array: Convert video to nd-array')
        probe = ffmpeg.probe(file)
        video_info = next(x for x in probe['streams'] if x['codec_type'] == 'video')
        self.width = int(video_info['width'])
        self.height = int(video_info['height'])

        out,err = (ffmpeg
                   .input(file)
                   .output('pipe:',format='rawvideo', pix_fmt='rgb24',**kwargs)
                   .run(capture_stdout=True))

        self.n_frames = int(len(out)/self.height/self.width/3)
        image_stack = np.frombuffer(out, np.uint8).reshape([-1, self.height, self.width, 3])
        return image_stack


    def crop_n_grayscale(self,video_array,
                         y=0,y_max=None,
                         x=0,x_max=None,
                         first_frame=None,
                         last_frame=None):
        print('---------------> Processing video:')
    
        ## Conditional for cropping video length
        if first_frame == None: first_frame = 0
        if last_frame == None: last_frame = video_array.shape[0]
        if y_max == None: y_max = video_array.shape[2]
        if x_max == None: x_max = video_array.shape[1]
        print('        - Cropping from frame %s to %s' % (first_frame,last_frame))
    
        ## Setting only frames and ROI to grayscale
        print('        - Converting to grayscale & cropping ROI to (%s x %s)' % (x_max-x,y_max-y))
        ch_1 = 0.2989 * video_array[first_frame:last_frame,y : y_max,x : x_max,0]
        ch_2 = 0.5870 * video_array[first_frame:last_frame,y : y_max,x : x_max,1]
        ch_3 = 0.1140 * video_array[first_frame:last_frame,y : y_max,x : x_max,2]
        clean_stack = ch_1.astype(float) + ch_2.astype(float) + ch_3.astype(float)

        print('  Final video array dimensions:',clean_stack.shape,'\n')
        return clean_stack
        
    def subtract_background(self,video_array=None,first_frame=0,last_frame=None):
        print('----> Calculate a background image')
        if last_frame == None: 
            last_frame = self.n_frames
#         if video_array == None: 
#             print('Exiting program. Needs an appropriate image stack to subtract background')
            
        background = np.median(video_array[first_frame:last_frame,:,:].astype(float), axis=0).astype(int)
        print('  subtract_background dimensions:     ', background.shape)
        print('')
        spot_stack = np.subtract(video_array,background)
        return spot_stack, background   


## Plots and views
    def view_ROI(self,image=None,border=True,x0=0,x1=None,y0=0,y1=None,**kwargs):
        print('---------------> view_ROI: View image with Region of Interest identified')

#         if image == None:
#             image = self.image_stack[0]

        plt.imshow(self.image_stack[0],cmap=cm.Greys_r, **kwargs)
        
        if border:
            if x1 == None: x1 = self.image_stack[0].shape[1]
            if y1 == None: y1 = self.image_stack[0].shape[0]
            plt.hlines(y0,x0,x1, color = 'r', alpha = .7)
            plt.hlines(y1,x0,x1, color = 'r', alpha = .7)
            plt.vlines(x0,y0,y1, color = 'r', alpha = .7)
            plt.vlines(x1,y0,y1, color = 'r', alpha = .7)
        return


    def display_images(self, cropped_converted, background,subtracted, frame=0, **kwargs):
        print('----> Displaying images')
        plt.figure(figsize = (6,8))
    
        ## Displaying the test frame image
        plt.subplot(311)
        plt.title('Cropped and converted, frame: %s' % str(frame))
        plt.imshow(cropped_converted[frame], cmap = cm.Greys_r, **kwargs)

        ## Displaying the background image
        plt.subplot(312)
        plt.title('Background image')
        plt.imshow(background, cmap = cm.Greys_r, **kwargs)

        ## Displaying the background subtracted, test frame image
        plt.subplot(313)
        plt.title('Subtracted background')
        plt.imshow(subtracted[frame], cmap = cm.Greys_r, **kwargs)

        plt.tight_layout()
        print('')
        return


    def image_metrics(self, spots, image, metric, colorbar=False, **kwargs):
        plt.title(metric)
        plt.imshow(image, cmap = cm.Greys_r)
        plt.scatter(spots.x,spots.y, c = spots[metric],
                    cmap = cm.coolwarm, **kwargs)
        if colorbar:
            plt.colorbar()
        plt.tight_layout()
        return


    def colored_hist(self, spots,metric,bins=None, predict_threshold=False, threshold=None):
        if bins == None: bins = 40
        try: threshold = int(threshold)
        except: pass
    
        set_cm = plt.cm.get_cmap('coolwarm')
        n, bin_assignments, patches = plt.hist(spots[metric],bins = bins)
        bin_centers = 0.5 * (bin_assignments[:-1] + bin_assignments[1:])
        col = bin_centers - min(bin_centers)
        col /= max(col)

        ## Plotting by color
        for c, p in zip(col, patches):
            plt.setp(p, 'facecolor', set_cm(c))
    
        ## Estimate auto-threshold
        if predict_threshold:
            _threshold = self.find_threshold(spots[metric],bins = bins)
            y_max = np.histogram(spots[metric],bins=bins)[0].max()
            plt.vlines(_threshold,0,y_max,color = 'gray',label='Auto')
        if isinstance(threshold, int) or isinstance(threshold, float):
            y_max = np.histogram(spots[metric],bins=bins)[0].max()
            plt.vlines(threshold,0,y_max,color = 'k', label='User-defined')
            if predict_threshold:
                plt.legend(frameon=False)
        return

    def spot_checker(self, spots,metrics=['signal'],image=None,**kwargs):
        ## Plotting each of the detector spot results over the image

        subplot = int(str(len(metrics)) + str(2) + str(0))
        if image==None: image = self.clean_stack[0]
        count = 0
        plt.figure(figsize=(4+image.shape[1]/150,len(metrics)*2))
        for i in range(0,len(metrics)):
            ## Defining spot metric and whether to auto-threshold
            col = metrics[i]
            if col=='signal': predict = True
            else: predict = False

            ## Drawing histogram, showing distribution of spots by metric
            count += 1
            plt.subplot(len(metrics),2,count)
            plt.title('Histogram for: %s' % col)
            self.colored_hist(spots,metric=col, bins = 40,predict_threshold=predict)
    
            ## Drawing image plot, colored by metric
            count += 1
            plt.subplot(len(metrics),2,count)
            plt.title('Spot overlay: %s' % col)
            self.image_metrics(spots,image, metric=col,**kwargs)

        plt.tight_layout()
        return
        
    def visualize_metrics(self,spots):
        metrics = ['ecc','mass','signal']
        self.spot_checker(spots,metrics=metrics, alpha=.1)
        return


## Spot finding
    def find_spots(self, stack, diameter = 3, quiet=True,**kwargs):
        ## Check diameter
        diameter = int(diameter)
        if ~diameter % 2: diameter = diameter + 1
    
        ## Option to silence output
        if quiet: tp.quiet()
    
        ## Detect spots
        spots = tp.batch(stack,diameter = diameter, **kwargs)

        ## Sorting DataFrame
        return spots[spots.raw_mass > 0].sort_values(by='frame')
        
        
    def particle_finder(self, invert=True, **kwargs): ## Break this into two functions, since parameter_testing function uses similar components
        print('---------------> particle_finder: Finds particles and formats resulting dataframe')
        df = self.find_spots(stack = self.spot_stack,
                                      quiet=True,invert=True,
                                      diameter=self.diameter,
                                      minmass=self.minmass,
                                      maxsize=self.maxsize)
        if df.shape[0] == 0:
            print('\n\nExiting program: No spots detected. Try modifying diameter, maxsize, or minmass parameters')
            raise SystemExit
        
        ## Round detector outputs
        df['x'] = df.x.round(2)
        df['y'] = df.y.round(2)
        df['t'] = [round(item/self.frame_rate,3) for item in df.frame]
        df['mass'] = df['mass'].astype(int)
        df['size'] = df.size.round(4)
        df['ecc'] = df.ecc.round(4)
        df['signal'] = df.signal.round(2)
        df['raw_mass'] = df['mass'].astype(int)
        df['ep'] = df.ep.round(1)
        df['True_particle'] = np.repeat(True,df.shape[0])#(df['signal'] > threshold)
        
        ## Saves DataFrame
        return df

    def find_threshold(self,x_array,bins=40):
        print('----> Selecting a threshold')
        x_array = np.histogram(x_array,bins = bins)[0]
        peaks = find_peaks(x_array)[0]
        prominences = peak_prominences(x_array,peaks)
        threshold = find_peaks(x_array, prominence=np.max(prominences))[0]
        print('    Threshold @',str(threshold[0]))
        print('')
        return threshold[0]

    def invert_Y(self,spots):
        return abs(spots.y - spots.y.max())
    
    def create_df_by_vial(self):
        '''
        Creates a dictionary with keys for vials and values for the DataFrame sliced by
        vial. It will also calculate the local linear regression for each vial and
        returns the DataFrame containing all the slopes and linear regression statistics
        for that vial.
        ----
        Inputs:
        df (DataFrame) : DataFrame sliced by vial
        vials (int) : Number of vials in video
        window (int) : Window size to calculate local linear regression
        vial_ID (str) : Vial-specific ID, taken from the first 'vvars' of naming convention
        ----
        Returns:
        result (dict) : DataFrame containing the local linear regression statistics
            by vial
        vial (dict) : Dictionary of DataFrames sliced by vial, keys are vials and values
            are DataFrames
        '''
        print('---------------> create_df_by_vial: Creates a DataFrame for each vial')
        self.vial,self.result = dict(),dict()


        for i in range(1,self.vials + 1):
            if i == 1: 
                self.vial[i] = self.df_filtered
            else: 
                self.vial[i] = self.df_filtered[self.df_filtered.vial==i]
            self.result[i] = self.local_linear_regression(self.vial[i]).iloc[0].tolist()
            self.result[i] = ['_'.join(self.vial_ID) + '_' + str(i)] + self.result[i]

            self.result[i][1:3] = [int(item) for item in self.result[i][1:3]]
            self.result[i][3:] = [round(item,4) for item in self.result[i][3:]]
        return


    def setting_final_slopes(self,df, pixel_to_cm, save=True):
        '''Taking most linear slope from each vial and saving it to a slopes.csv file'''
        print('---------------> setting_final_slopes: ')
        cols = ['vial_ID','first_frame','last_frame','slope','intercept','r_value','p_value','std_err']
        self.df_slopes = pd.DataFrame.from_dict(df,orient='index',columns = cols)
        if self.convert_to_cm:
            self.df_slopes['slope'] = self.df_slopes['slope'] / int(pixel_to_cm)
            self.df_slopes['slope'] = self.df_slopes['slope'].round(4)

        if save:
            self.df_slopes.to_csv(self.name_nosuffix + '.slopes.csv',index=False)
        return self.df_slopes


    def bin_vials(self,df, vials):
        '''Bin spots into vials. Function takes into account all points along the x-axis, 
        and divides them into specified number of bins based on the min and max
        points in the array.
        ----
        Note: Number of vials is specified in <file>.cfg'''
        print("---------------> bin_vials: assign spots to a vial")
        bin_lines = pd.cut(df.x,vials,include_lowest=True,retbins=True)[1]
        spot_assignments = pd.cut(df.x,bins=bin_lines,labels=range(1,vials+1))
        ## Important if a middle vial disappears
        counts = np.unique(spot_assignments,return_counts = True)
        for v,c in zip(counts[0],counts[1]):
            if c == 0:
                print('Warning: vial',v,'is empty and cannot be evaluated')
        return bin_lines, spot_assignments


    def local_linear_regression(self, df, method = 'max_r'):
        '''Performs a local linear regression, using a user-defined sliding window'''
        print("---------------> local_linear_regression")
        result_list,result = [],pd.DataFrame()
        llr_columns = ['first_frame','last_frame','slope','intercept','r','pval','err']

        for i in range(len(self.image_stack) - self.window):
            start,stop = int(i),int(i+self.window)
            df_window  = df[(df.frame >= start) & (df.frame <= stop)]

            ## Testing if there are enough frames in the slice
#             if len(df_window.frame.unique()) < self.window-2: # Kept getting errors checking to see if there were frames with no spots. Going with groupby function instead
            if df_window.groupby('frame').y.count().min() == 0:
                print('Issue with number of frames with flies detected vs. window size')
                print(i,i+self.window,len(df_window.frame.unique()))
                continue
                        
            ## Performs local linear regression
            try: 
                ## Grouping points by frame
                _frame = df_window.groupby('frame').frame.mean()
                _pos  = df_window.groupby('frame').y.mean()

                ## Performing linear regression on subset and formatting output
                _result = linregress(_frame,_pos)
                _result = [start,stop] + np.hstack(_result).tolist()

                ## If slope is not significantly different from 0, then set slope = 0
                if _result[-2] >= 0.05: _result[2] = 0
            except: _result = [start,stop] + [np.nan,np.nan,np.nan,np.nan,np.nan]

            result_list.append(_result)
        
        result = pd.DataFrame(data=result_list,columns=llr_columns)

        if method == 'max_r':
            result = result[result.r == result.r.max()]
        elif method == 'min_err':
            result = result[result.err == result.err.min()]
        else:
            print("Unrecognized method, chose either 'max_r' to select the window with the greatest R or 'min_err' to select the window with the lowest error")
        return result
        
        
    def step_1(self):
        '''Words'''
        print('-- [ Step 1 ] Cleaning and format image stack')
        y = self.y
        y_max = int(y + self.h)
        x = self.x
        x_max = int(x + self.w)
        stack = self.image_stack
        self.clean_stack = self.crop_n_grayscale(stack,
                         y=y, y_max=y_max,
                         x=x, x_max=x_max)
        print(self.clean_stack.shape)

        self.spot_stack,self.background = self.subtract_background(video_array=self.clean_stack)
#         plt.figure()
#         self.display_images(self.clean_stack,self.background,self.spot_stack,frame=20)
#         plt.savefig(self.name_nosuffix + '.processed.png', dpi=100)
#         plt.close()
        return


    def step_2(self):
        print('-- [ Step 2a ] Detecting flies')
        self.df_big = self.particle_finder(minmass=self.minmass,diameter=self.diameter,
                                            maxsize=self.maxsize, invert=True)
        self.df_big.to_csv(self.name_nosuffix + '.raw.csv', index=None)
        print('                -->',self.path_data.split('/')[-1])

        ## Auto-detecting threshold
        if self.threshold == 'auto':
            self.threshold = self.find_threshold(self.df_big.signal)
        
        self.df_big['True_particle'] = self.df_big.signal.map(lambda x: x >= self.threshold)
        if self.df_big[self.df_big.True_particle].shape[0] == 0:
            print('\n\nExiting program: Threshold limit is too restrictive, no spots left after filtering')
            raise SystemExit
        
        print('-- [ Step 2b ] Detecting flies, filtering and binning into vials')        
        ## Filtering spots and pruning unnecessary columns
        self.df_filtered = self.df_big[self.df_big.True_particle]
        self.df_filtered = self.df_filtered.drop(['ecc','signal','ep','raw_mass','mass','size','True_particle'],axis=1)
    
        ## Binning spots into vials
        self.bin_lines, self.df_filtered['vial'] = self.bin_vials(self.df_filtered,vials = self.vials)
        self.df_filtered = self.df_filtered.sort_values(by=['vial','frame','y','x'])

        print('-- [ Step 2c ] Adding experimental variables to DataFrame\n')
        for item in self.file_details.keys():
            self.df_filtered[item] = np.repeat(self.file_details[item],self.df_filtered.shape[0])
        
        return

    def step_3(self):
        print('-- [ Step 3 ] Visualize spot metrics')

        plt.figure()
        self.view_ROI(border = True,
                        x0 = self.x, x1 = self.x + self.w,
                        y0 = self.y, y1 = self.y + self.h)
        plt.savefig(self.name_nosuffix + '.ROI.png',dpi=100)
        plt.close()
                    
        self.visualize_metrics(self.df_big)
        plt.savefig(self.name_nosuffix + '.spot_check.png',dpi=300)
        plt.close()
        return
        
    def step_4(self):
        print('-- [ Step 4 ] Calculating local linear regressions')
        
        ## Invert the y-axis
        self.df_filtered['y'] = abs(self.df_filtered.y.max() - self.df_filtered.y)
        self.df_filtered['y'] = self.df_filtered.y.round(2)

        ## Save the filtered DataFrame
        self.df_filtered.to_csv(self.name_nosuffix+'.filtered.csv', index=False)
        print('                -->',self.path_filter.split('/')[-1])
        return
        
    def step_6(self,gui=False):
        print('-- [ Step 6 ] Plotting data')
        plt.figure(figsize=(10,8))
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
#         fig.suptitle(self.name)

        ## Finding the frames that flank the most linear portion 
        ##    of the y vs. t curve for all points, not just by vials
        begin = self.local_linear_regression(self.df_filtered).iloc[0].first_frame.astype(int)
        end = begin + self.window

        print("Filtering")
        ## Filtering for true points, then binning into vials
        spots = self.df_big[self.df_big.True_particle]
        self.bin_lines, self.df_big['vial'] = self.bin_vials(self.df_big, self.vials)

        ## Plot: Image of 'begin' frame with spots
        print("Plot 1")
        self.image_plot(df = spots,frame = begin, ax=ax1) # ADAM...issue if only one vial
        ## Plot: Image of 'end' frame with spots
        print("Plot 2")
        self.image_plot(df = spots, ax=ax2, frame = int(str(end)))
        ## Plot: Image of 0th frame with all spots from video
        print("Plot 3")
        self.image_plot(df = spots,ax=ax4, frame = None)
        print("Create DF")                  
        ## Breaking DataFrame into vials
        self.create_df_by_vial()
        print("Plot 4")
        self.loclin_plot(ax=ax3)
        fig.tight_layout()
        plt.savefig(self.name_nosuffix + '.diagnostic.png',dpi=500, transparent = True)
        plt.close()
        print('                --> Saved:',self.path_plot.split('/')[-1])
        
        if gui:
            return ax
        else:
            return
    
    def step_7(self):
        print('-- [ Step 7 ] Setting up slopes file')
        slope_columns = ['vial_ID','first_frame','last_frame','slope','intercept','r_value','p_value','std_err']
        self.df_slopes = pd.DataFrame.from_dict(self.result,orient='index',columns = slope_columns)
        if self.convert_to_cm:
            self.df_slopes['slope'] = self.df_slopes.slope.transform(lambda x: x * (self.pixel_to_cm / self.frame_rate)).round(6)
        
        for item in self.file_details.keys():
            self.df_slopes[item] = np.repeat(self.file_details[item],self.df_slopes.shape[0])
    
        ## Saving slope file
        self.df_slopes.to_csv(self.name_nosuffix + '.slopes.csv',index=False)
        print('                --> Saved: %s \n' % self.path_slope.split('/')[-1])
        return
        
    def image_plot(self,df,frame=None,ax=None,ylim=[0,1000]):
        ## Can probably replace this with the image plot function in the test_parameters function below
        ## Assigning plotting parameters based on frame
        try: frame = int(frame)
        except: frame = None

        ## Determine which frame to plot
        if frame == None:
            frame = 0
            alpha = 0.01
            ax.set_title('All x,y-points throughout video')
        elif type(frame) == int:
            df = df[(df.frame == frame)]
            alpha = .25
            ax.set_title("Frame: %s" % frame)
        else:
            print("Chose a frame value in integer form, or 'None'")

        ## Issue with plotting if last frame in stack
        if frame == self.n_frames:
            image = self.clean_stack[frame-1]
        else:
            image = self.clean_stack[frame]

        ## Plotting image
        ax.imshow(image,cmap=cm.Greys_r,origin='upper')
        
        ## Plotting vertical bin lines
        ax.vlines(self.bin_lines,0,image.shape[0],alpha = .3)  

        ## Plotting point overlay
        df['vial'] = self.bin_vials(df,self.vials)[1] ## Issue if only one vial?
        df = df.sort_values(by='vial')
        if self.vials > 1:
            ax.scatter(df.x,df.y,s=30,alpha=alpha,
                    c = df['vial'], cmap=self.discrete_cmap(vials = self.vials))
        else:
            ax.scatter(df.x,df.y,s=30,alpha=alpha,
                    color = self.color_list[0])
       
        ## Getting rid of axis labels
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        return ax
    
    def loclin_plot(self,ax = None):
        print("---------------> loclin_plot")
        def two_plot(df,color,label,first, last, ax=None):
            print("---------------> two_plot")
            
            x = df.groupby('frame').frame.mean()
            y = df.groupby('frame').y.mean()
            ax.plot(x,y, alpha = .35, color = color,label='') 

            df = df[(df.frame >= first) & (df.frame <= last)]
            x = df.groupby('frame').frame.mean()
            y = df.groupby('frame').y.mean()
            ax.plot(x,y,color = color,label = label)
            return
        
        if len(self.result) > 1:
            for V in range(1,self.vials+1):
                l = 'Vial '+str(V)
                c = self.color_list[V-1]
                _details = self.result[V]
                frames = _details[1],_details[2]
                two_plot(self.vial[V],color = c,label = l,
                         first = _details[1],last = _details[2],ax=ax)
        else:
            for V in range(1,self.vials+1):
                l = 'Vial 1'
                c = self.color_list[V-1]
                _details = self.result[V]
                frames = _details[1],_details[2]
                two_plot(self.vial[V],color = c,label = l,
                         first = _details[1],last = _details[2],ax=ax)

        ax.legend(loc=2, frameon=False, fontsize='x-small')

        label = 'pixels'
        if self.convert_to_cm: label = 'cm'
        labels = ['Cohort climbing kinematics', 'Frames','Mean y-position (%s)' % label]
        ax.set(title=labels[0],xlabel = labels[1],ylabel=labels[2])
        return

    def discrete_cmap(self, vials = None):
        ## Modified from: https://gist.github.com/jakevdp/91077b0cae40f8f8244a
        '''Creates a list of colors from a color map, given a specified number of vials'''
        if vials == None:
            vials = 1
            base_cmap = 'tab10'
        elif self.vials <= 10: base_cmap = 'jet'
        elif self.vials <= 20: base_cmap = 'tab20'
        elif self.vials > 20: base_cmap = 'jet'

        base = plt.cm.get_cmap(base_cmap)
        
        ## If there is an odd number of vials, adjust for having a too-light 'middle' color
        if vials%2 and vials > 1:
            self.color_list = base(np.linspace(0, 1, vials+1))
            row_to_drop = len(self.color_list)//2
            self.color_list = np.delete(self.color_list, row_to_drop,axis = 0)
        else:
            self.color_list = base(np.linspace(0, 1, vials))
        
        self.cmap_name = base.name + str(self.vials)
        return LinearSegmentedColormap.from_list(self.cmap_name, self.color_list, vials)

    
    #####    #####    #####    #####    #####    #####    #####    #####    #####
    ## Commands used only in GUI
    #####    #####    #####    #####    #####    #####    #####    #####    #####
    def parameter_testing(self, variables, axes):
        print('-------> parameter_testing')
        self.load_for_gui(variables)
        self.step_1() # Load video
        self.step_2() # Crop video
        self.step_3() # Visualize other spots
        
        spots = self.df_big
        SPOTS = spots[spots.True_particle]
        bin_lines,SPOTS['vial'] = self.bin_vials(SPOTS,self.vials)
        self.discrete_cmap(vials=self.vials)
        SPOTS.loc[:,'color'] = SPOTS.vial.map(dict(zip(range(1,self.vials+1),self.color_list)))
        print(SPOTS.head())
        frame_0 = self.frame_0
        frame_n = self.frame_n
        bins=20
        
        ## Setting plot for background image
        axes[0].set_title("Background Image")
        axes[0].imshow(self.background,cmap=cm.Greys_r)
        axes[0].scatter([0,self.w],[0,self.h],alpha=0,marker='.')

        ## Setting plots for scatterplot overlay
        frame=40
        axes[1].set_title('Frame:'+str(frame_0))
        axes[1].imshow(self.clean_stack[frame_0], cmap = cm.Greys_r)
        axes[1].scatter(spots[(spots.frame==frame_0) & (~spots.True_particle)].x,
                        spots[(spots.frame==frame_0) & (~spots.True_particle)].y, 
                        color='w',marker = 'x', alpha = .5)
        a = axes[1].scatter(SPOTS[SPOTS.frame==frame_0].x,
                            SPOTS[SPOTS.frame==frame_0].y, 
                            color = SPOTS[SPOTS.frame==frame_0].color,
                            marker ='o',alpha = .8)
        a.set_facecolor('none')
        axes[1].vlines(bin_lines,0,spots.y.max(),color='w')
        
        ## Fly counts
        axes[5].plot(SPOTS.groupby('frame').frame.mean(),
                     SPOTS.groupby('frame').frame.count(), 
                     label = 'Fly count', color = 'g',alpha = .5)
        axes[5].hlines(np.median(SPOTS.groupby('frame').frame.count()),
                       spots.frame.min(),spots.frame.max(),
                       linestyle = '--',alpha = .5, 
                       color = 'gray', label = 'Median no. flies')
        axes[5].set(title = 'Flies per frame',
                        ylabel='Flies past threshold',
                        xlabel='Frame') 
        axes[5].legend(frameon=False, fontsize = 'small')

        ## Mass histogram
        axes[3].set_title('Mass Distribution')
        axes[3].hist(spots.mass,bins = bins)
        y_max = np.histogram(spots.mass,bins=bins)[0].max()
        axes[3].vlines(self.threshold,0,y_max)

        ## Signal histogram
        axes[4].set_title('Signal Distribution')
        axes[4].hist(spots.signal,bins = bins)
        y_max = np.histogram(spots.signal,bins=bins)[0].max()
        axes[4].vlines(self.threshold,0,y_max)

        ## Setting plots for local linear regression
        df = SPOTS.sort_values(by='frame')
        df['y'] = self.invert_Y(df)
        for V in range(1,self.vials + 1):
            label = 'Vial '+str(V)
            color = self.color_list[V-1]
            _df = df[df.vial == V]
            
            begin = self.local_linear_regression(_df).iloc[0].first_frame.astype(int)
            end = begin + self.window        
            frames = begin,end
            
            axes[2].plot(_df.groupby('frame').frame.mean(),
               _df.groupby('frame').y.mean(),alpha = .35, color = color,label='') 
            _df = _df[(_df.frame >= begin) & (_df.frame <= end)]

            axes[2].plot(_df.groupby('frame').frame.mean(),
               _df.groupby('frame').y.mean(),color = color, label = label)

        labels = ['Mean vertical position over time','Mean vertical position','Frame']
        axes[2].set(title = labels[0], ylabel=labels[1],xlabel=labels[2]) 
        axes[2].legend(frameon=False, fontsize='x-small')            
        return
#         