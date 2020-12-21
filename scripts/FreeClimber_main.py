#!/usr/bin/env python
# -*- coding: utf-8 -*- 

## File name : FreeClimber_main.py
## Created by: Adam N. Spierer
## Date      : May 2020
## Purpose   : Command line interface wrapper for FreeClimber


## Version number
version='0.4.0'
doi =  'https://doi.org/10.1242/jeb.229377' ## Link to published paper

## Importing external package(s)
import os
import sys
import argparse
from time import time,ctime
from datetime import datetime
from pandas import read_csv, concat
from numpy import sort, repeat, unique
import matplotlib.pyplot as plt
from matplotlib.cm import Greys_r

## Importing local module(s)
import detector as detector

class FreeClimber(object):
    def __init__(self, config_file):
        '''Initializing the detector with the configuration file and arguments'''
        self.config_file = config_file
        self.args = define_argument_parser()
        if self.args.debug: print('FreeClimber.__init__')
        
        ## Load main parameters and get the list of files to process
        self.load_parameters()
        self.get_filelist()
        ## Insert 1 -- variable check

        ## Basic variables
        self.count = 0
        self.first_run = True
        return
            
    def load_parameters(self):
        '''Loads parameters for the FreeClimber object'''
        if self.args.debug: print('FreeClimber.load_parameters')
        
        ## Read in parameters from configuration file
        var_list = []
        with open(self.config_file,'r') as f:
            variables = f.readlines()#.split('\n')
        f.close()

        ## Filter lines with '#', ' ', and carriage returns
        variables = [item.rstrip() for item in variables if not item.startswith(('#',' ','\n'))]
        
        ## Assign variables to the detector, pass if unable to.
        for item in variables:
            try: exec('self.'+item)
            except: pass
        return

    ## Reading file with video paths for --process_custom argument
    def read_custom(self, file):
        '''Reads in process file (.prc) containing individual video file paths to process.
        File path and type were confirmed earlier in get_filelist.
        ----
        Inputs:
          file (str): Path to custom file (.prc)
        ----
        Returns:
          _lines (list): List of valid file paths read from custom file (.prc)
        '''
        if self.args.debug: print('FreeClimber.read_custom')
        
        ## Read in lines from custom file
        with open(file) as f:
            lines = [line.rstrip() for line in f]
        f.close()

        ## Checks individual files are real
        _lines = [vid for vid in lines if os.path.isfile(vid)]
        print('----> %s of %s specified files have valid paths' % (len(lines),len(lines)))
        return _lines
    
    def get_filelist(self):
        '''Gets the list of files to process. Defaults to all'''
        if self.args.debug: print('FreeClimber.get_filelist')
        
        ## Checks for undone argument
        if self.args.process_undone:
            self.file_list = self.file_walker(self.path_project, endswith=self.file_suffix, undone = True)
            return
        
        ## Checks for custom files, with specific error messages w/ invalid file path or type
        elif self.args.process_custom != False:
            if os.path.isfile(self.args.process_custom):
                if self.args.process_custom.endswith('.prc'):
                    self.file_list = self.read_custom(file = self.args.process_custom)
                    return
                else:
                    print('!Not valid!', self.args.process_custom,'\n')
                    print("--> Exiting program: Custom file lacks '.prc' suffix")
                    raise SystemExit
            else:
                print('!Not valid!', self.args.process_custom,'\n')
                print("--> Exiting program: Custom file path is invalid")
                raise SystemExit           
                         
        ## Defaults to processing all files
        else:
            self.file_list = self.file_walker(self.path_project, endswith=self.file_suffix)

    def print_new_project(self):
        '''Prints out a header for the new project'''
        if self.args.debug: print('FreeClimber.print_new_project')
        
        print("##\n## Analyzing videos in: %s \n##" % self.path_project )
        print("## Processing %s videos\n##\n##" % (str(len(self.file_list))))
        return

    def file_walker(self,folder = None, endswith = None, undone = False):
        '''From a specified parent folder, find all child files with the specified suffix
        ----
        Inputs:
          folder (str): Parent folder to search through
          endswith (str): Suffix of a common file type
          undone (bool): False finds all files with the 'endswith' suffix. 
                         True does the same, but excludes '.slopes.csv' files (processed)
        ----
        Returns:
          _list (list): sorted list of all file paths with a common suffix in a parent folder
        '''
        if self.args.debug: print('FreeClimber.file_walker')
        
        _list1,_list2 = [],[]
        for root, dirs, files in os.walk(folder):
            for name in files:
                if name.endswith(endswith):
                    _list1.append((os.path.join(root, name)))
                if undone:
                    if name.endswith('.slopes.csv'):
                        _list2.append(os.path.join(root, name[:-11])+'.'+endswith)

        ## Return a sorted list of all files with the file suffix
        if undone == False:
            _list = sorted(unique(_list1))
            return _list

        ## Return a sorted list of all undone files            
        if undone:
            _list1,_list2 = set(_list1),set(_list2)
            _list = list(_list1.difference(_list2))
            _list = sorted(	_list)
            if len(_list) == 0:
                print('All files previously processed, re-evaluate your inputs if this message is a surprise.')
            return 	_list

    def timer(self, time_begin):
        '''Timer for measuring each video's processing time'''
        if self.args.debug: print('FreeClimber.timer')
        
        time_elapsed = time()-time_begin
        time_elapsed = str(int(time_elapsed//60))+" min : " + str(int(time_elapsed%60))+" seconds"
        print("Time elapsed: ",time_elapsed)
        print('='*72)
        return
        
    def print_new_video(self,name):
        '''Prints header for a new video with file name, how many processed, and remainder'''
        if self.args.debug: print('FreeClimber.print_new_video')

        print('='*72)
        line_to_print = '== [%s || %s] ' % (self.count,len(self.file_list)), self.name
        line_to_print = ('').join(line_to_print)

        if len(line_to_print) <= 80: print(('').join(line_to_print), '=' * (71 - len(line_to_print)))
        else: print(('').join(line_to_print))
        print('='*72)
        return

    def process(self,video_file,variables, config_file = None, gui=False):
        '''Executes the steps in the detector object'''
        if self.args.debug: print('FreeClimber.process')

        if config_file == None:
            print('Requires a valid configuration file')
        else:
            self.print_new_video(video_file)
            d = detector.detector(video_file = video_file, config_file = config_file, debug = self.args.debug)
            d.step_1(gui = self.args.optimization_plots) # Crops and formats the video
            d.step_2() # DataFrame creation and manipulation (df_big and df_filtered)

            d.step_3(gui = self.args.optimization_plots) # Visualizes spot metrics
            d.step_4()# Filters and processes data detected points

            d.step_5() # Calculates local linear regressions
            d.step_6(gui = self.args.optimization_plots) # Creating diagnostic/other plots 
            d.step_7() # Writing the video's slope file
            self.first_run = False
        return

    def concat_slopes(self):
        '''Concatenate the .slopes.csv files into a single results.csv file in path_project folder'''
        if self.args.debug: print('FreeClimber.concat_slopes')        

        ## Finds all files with the .slopes.csv suffix
        print("\nFinal step: Concatenating slope files")
        slope_files = self.file_walker(self.path_project,endswith='.slopes.csv')
        to_concat = [slope_file for slope_file in slope_files]
        if self.args.debug: print(to_concat)
        
        ## Concatenates contents of files into a single DataFrame
        print("    - Concatenating",len(to_concat),"files")
        self.slopes = concat([read_csv(item) for item in to_concat])
        
        ## Saves results.csv file to the path_project folder (specified in the configuration file)
        self.path_result = self.path_project+'results.csv'
        print("    - Saving:", self.path_result)
        self.slopes.to_csv(self.path_result,index=False)
        return

    def create_log_header(self):
        '''Creates a header for the completed and skipped log files'''
        if self.args.debug: print('FreeClimber.create_log_header')
    
        try:
            os.mkdir(self.path_project + 'log/')
        except:
            pass
        self.path_completed = self.path_project + 'log/completed.log'
        self.path_skipped = self.path_project + 'log/skipped.log'
        path_list,text_list = [self.path_completed,self.path_skipped],['completed','skipped']
        time_stamp = str(ctime())
                
        ## For future release
#         if self.review_R != None:
#             self.path_review = self.path_project + '/review_at_R_%s/review.log' % self.review_R
#             path_list.append(self.path_review)
#             text_list.append('review')
    
        for path,text in zip(path_list,text_list):
            print('Creating',text,'.log file: ',path)
            with open(path, 'w') as f:
                print('## FreeClimber ##', file=f)
                print('##',file=f)
                print('## Generated from configuration file: '+self.config_file, file = f)
                print('## Run on ' + time_stamp, file = f)
                print('##',file = f)
                print('## Files %s:' % text, file = f)

## For future release
#                 if self.review_R != None: print('## Files %s:' % text, file = f)
#                 else: print('## Files for review @ R = %s' % self.review_R, file = f)

            f.close()
        return

    def log_video(self, file_name, completed=True): #, review=None):
        '''Records which videos were processed and which were skipped.
        ----
        Inputs:
          file_name (str): file path of video being skipped
          completed (bool): True if successfully processed, False if skipped        
        ----
        Returns:
          None'''
        if self.args.debug: print('FreeClimber.log_video :: completed =',completed)
    
        ## Naming a variable
        if completed: path = 'log/completed.log'
        else:  path = 'log/skipped.log'

## For future release
#         elif review != None: path = '/review_at_R_%s/review.log' % self.review_R
        
        ## Appends the video_file path to the log file
        if ~completed:
            print('Appending to',path + ': ',file_name)
        with open(self.path_project + '/' + path,'a') as f:
            print(file_name,file = f)
        f.close()
        return
        
    def print_closing(self):
        '''Print statements when program is complete. Need to add in a counter to display 
        how many were processed'''
        if self.args.debug: print('FreeClimber.print_closing')        
        
        print('\nVideo processing complete!\n')
        print('Total videos: %s' % len(self.file_list))
## Upcoming added functionality to output the number of videos processed vs initially input
#         print('Processed   : %s')
#         print('Unprocessed : %s ... see ## FILE PATH ##')
        return


def define_argument_parser():
    '''Defines arguments to be parsed, via argparse module.
    
    The '--file' flag is required and points to a configuration file (.cfg)
    created from the GUI. 
    
    Files to process are specified for all (--process-all), undone (--process_undone), 
    or a custom list (--process_custom) of videos. If a custom list, include the path to 
    a process file (.prc) containing all files to process, each on a separate line. The
    gather_files.py script can help generate one, or the log/skipped.log file will
    contain all paths that were passed over due to technical difficulties.
    
    The '--no_concat' flag will prevent the final concatenation step, which is useful if
    users want to analyze individual videos but not overwrite the results.csv file.
    
    The '--optimization_plots' flag will create the optimization plots generated by the 
    GUI. These include files with suffixes: ROI.png, spot_check.png, and processed.png.
    
    ## For future release
    The '--review_R' flag and argument will create a list of files with vials that have
    a regression coefficient (R) value that is less than a predefined threshold 
    (default = None) and save a lower resolution copy of the diagnostic plot for review. 

    ----
    Inputs:
      None
    ----
    Returns:
      args (object): Namespace object containing the flags and arguments passed to program
    '''    
    ## Set up the command line parser
    parser = argparse.ArgumentParser(prog='FreeClimber',
                                    description='Calculates the climbing velocity of flies in a negative geotaxis assay',
                                    #  usage='%(prog)s [options] path',
                                    epilog='For documentation and a tutorial, see https://github.com/adamspierer/FreeClimber',
                                    allow_abbrev=False)

    ## Specifying configuration file
    parser.add_argument('--config_file', 
                        type=str, 
                        required=True, 
                        help="Path to configuration file (ends with '.cfg')")

    ## Batch processing: specifying which files to process
    group_method = parser.add_mutually_exclusive_group(required=False)
    group_method.add_argument('--process_all', 
                        default=False, 
                        action='store_true', 
                        help="Process all files")

    group_method.add_argument('--process_undone', 
                        default=False,
                        action='store_true', 
                        help="Process previously unprocessed files")
    group_method.add_argument('--process_custom', 
                        default=False, 
                        type=str,
                        help="Process custom list of files, must include file path as argument")

    ## Batch processing: additional arguments
    parser.add_argument('--no_concat', 
                        required=False, 
                        default=False, 
                        action='store_true',
                        help="Use this flag to prevent results files from concatenating")

    ## Generate all possible plots
    parser.add_argument('--optimization_plots', 
                        required=False, 
                        default=False, 
                        action='store_true',
                        help="Creates the detector optimization plots with spot metrics for each video")

    ## For future release
    ## Specify regression coefficient for secondary review
#     parser.add_argument('--review_R', 
#                         required=False, 
#                         default=None,
#                         type=float,
#                         help="Regression coefficient threshold for flagging noiser velocity curves")
    
    ## Debug printing
    parser.add_argument('--debug', 
                        required=False, 
                        default=False, 
                        action='store_true',
                        help="Prints each function's name as it is run--not very pretty")

    ## Final parser specifications
    parser.version = '1.0' ## Program version
    parser.add_argument('-v','--version', 
                        action='version')
    args = parser.parse_args()
    if args.debug: print(args)
    return args


## Checking argparse arguments
def check_config(args):
    '''Checks arguments passed to parser and kills program if invalid path.
    ----
    Inputs:
      args (list): Arguments passed to parsed
    ----
    Returns:
      file (str): File path if it is entered and valid, or None if not.
    '''
    if args.debug: print('__main__.check_config')

    ## Confirm file path is valid, exit if not
    if os.path.isfile(args.config_file): pass
    else:
        print('!Not valid!', args.config_file,'\n')
        print("--> Exiting program: Invalid file path, not a configuration file, or file lacks '.cfg' suffix")
        raise SystemExit
    return args.config_file
    
    
def startup():
    '''Function prints top line and runs argument parsing function.'''
    line_length = 72
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def print_line(line,line_length):
        '''Formats line to print'''
        if len(line) <= line_length: string = line + '#'*(line_length-len(line))
        else: string = line
        print(string)
        return

    ## Lines to print
    line0 = '#'*line_length
    line1 = '## FreeClimber v.%s ' % str(version)
    line2 = '## Please cite: %s' % doi
    line3 = "## Beginning program @ %s " % str(now)
    line4 = line0
    
    ## Printing formated lines
    print('\n')
    for item in range(5):
        print_line(eval('line'+str(item)),line_length)
    return         

def main():
    '''Main loop for FreeClimber'''
    ## Printing start sequence
    startup()
    
    ## Defining argparse commands
    args = define_argument_parser()

    ## Checking input for valid input file(s)
    config_file = check_config(args)

    ## Set up FreeClimber object and load parameters
    fc = FreeClimber(config_file = config_file)
    fc.create_log_header()

    for File in fc.file_list:
        if args.debug:
            print(File)
            if 1==1:
                t0 = time()
                fc.count += 1
                fc.name = os.path.split(File)[-1]
                fc.process(video_file = File,variables = None, config_file = fc.config_file)
                fc.timer(t0)
                fc.log_video(completed=True, file_name = File)
            else:
                fc.log_video(completed=False, file_name = File)    
            
        else:
            try:
                t0 = time()
                fc.count += 1
                fc.name = os.path.split(File)[-1]
                fc.process(video_file = File,variables = None, config_file = fc.config_file)
                fc.timer(t0)
                fc.log_video(completed=True, file_name = File)
            except:
                fc.log_video(completed=False, file_name = File)    

    ## Concatenate slopes of all .slopes.csv files into a single, results.csv file
    if args.no_concat == False:
        fc.concat_slopes()
        
    fc.print_closing()
    return
        
if __name__ == '__main__':
    main()
    exit()


################################################################
## Additional functionality to add in for later releases
    
## Insert 1
#    fc.check_variables() # Not ready for this release
#
#     def check_variables(self):
#     
#         ## Making sure diameter is an odd number
#         if int(self.diameter) % 2 == False:
#             self.diameter = int(self.diameter) + 1
#             
#         ## Check other variables...
#         return
