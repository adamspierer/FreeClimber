#!/usr/bin/env python
# -*- coding: utf-8 -*- 

## Upcoming features:
## [ ] Argparse options for processing {all | undone | custom} files
## [ ] Test statements to confirm variables are entered properly
## [ ] Log files for all outputs
## [ ] Log files for videos skipped

## Version number
version='0.3'

## Importing external packages
import os
import detector as detector
import sys
import argparse
from time import time
from datetime import datetime
from numpy import sort, repeat
from matplotlib.cm import Greys_r
from numpy import unique
import pandas as pd
import matplotlib.pyplot as plt


def load_variables(config_file):
    ## Finding the path_project
    if os.path.isfile(config_file):
        pass
    else:
        print('\n\nExiting program. Invalid configuration file')
        raise SystemExit

    with open(config_file,'r') as f:
        print('## Reading in configuration file:', config_file,'\n')
        variables = f.read().split('\n')
    f.close()

    ## Importing configuration file variables
    variables = [item for item in variables if ~item.startswith(' ') or item.startswith('#')]
    return variables


class FreeClimber(object):
    def __init__(self, cfg_file):
        self.cfg_file = cfg_file
        self.args = define_argument_parser()
    
            
    def load_parameters(self):
        var_list = []
        with open(self.cfg_file,'r') as f:
            variables = f.read().split('\n')
        f.close()
        variables = [item for item in variables if not item.startswith('#') or not item.startswith(' ') or not item.startswith('\n')]
        for item in variables:
            try:
                exec('self.'+item)
            except:
                pass
        return
    
    ## Insert 4

    def file_walker(self,folder = None,endswith=None):
        '''
        From a specified parent folder, find all child files with the specified suffix
        '''
        _lst=[]
        for root, dirs, files in os.walk(folder):
            for name in files:
                if name.endswith(endswith):
                    _lst.append((os.path.join(root, name)))
        return sorted(unique(_lst))

    def undone(self,folder=None):
        '''
        Used in conjunction with file_walker, this function finds all the files that have
        yet to be processed.
        '''
        _lst1,_lst2 = [],[]
        for root,dirs,files in os.walk(folder):
            for name in files:
                if name.endswith(self.file_suffix):
                    _lst1.append(os.path.join(root, name[:-5]))
                if name.endswith('slopes.csv'):
                    _lst2.append(os.path.join(root, name[:-4]))

        _lst1,_lst2 = set(_lst1),set(_lst2)

        todo = list(_lst1.difference(_lst2))
        todo = [item+'.'+self.file_suffix for item in todo]
        return sorted(todo)

    def timer(self, time_begin):
        time_elapsed = time()-time_begin
        time_elapsed = str(int(time_elapsed//60))+" min : " + str(int(time_elapsed%60))+" seconds"
        print("Time elapsed: ",time_elapsed)
        print('='*72)
        return

    def setup(self):
        '''Basic variables to begin'''
        self.count = 0
        self.first_run = True
        return
        
    def print_new_video(self,name):
        print('='*72)
        line_to_print = '== [%s || %s] ' % (self.count,len(self.file_list)), self.name
        line_to_print = ('').join(line_to_print)
        print(('').join(line_to_print), '=' * (71 - len(line_to_print)))
        print('='*72)
        return

    def print_new_project(self):
        print("##\n## Analyzing videos in: %s \n##" % self.path_project )
        print("## Processing %s videos\n##\n##" % (str(len(self.file_list))))
        return

    def process(self,File,variables, gui=False):
        self.print_new_video(self.name)
        d = detector.detector(video_file = File, config_file = self.cfg_file)
        d.step_1()
        d.step_2()
        d.step_3(gui=False) # Creating spot check plots for gui
        d.step_4() # Calculating local linear regression
        d.step_5() # Plotting data
        d.step_6() # Setting up slopes file
        
        self.first_run = False
        return

    def concat_slopes(self):
        print("\nFinal step: Concatenating slope files")
        slope_files = self.file_walker(self.path_project,endswith='.slopes.csv')
        to_concat = [slope_file for slope_file in slope_files]
        print("    - Concatenating",len(to_concat),"files")
        self.slopes = pd.concat([pd.read_csv(item) for item in to_concat])
        self.path_result = self.path_project+'/results.csv'
        print("    - Saving:", self.path_result)
        self.slopes.to_csv(self.path_result,index=False)
        return
    
    ## Insert 6
    
    def print_closing(self): ## ADAM: FINISH THIS
        '''Print statements when program is complete. Need to add in a counter to display how many were processed'''
        print('\nVideo processing complete!\n')
#         print('Total videos: %s' % )
#         print('Processed   : %s')
#         print('Unprocessed : %s ... see ## FILE PATH ##')
        return


def define_argument_parser():
    '''Defines arguments to be parsed, via argparse module. 
    Given a method argument, 'gui' or 'batch', different downstream arguments are required.
    
    The '--file' flag is required and points to a configuration file (.cfg)
    created from the GUI. 
    
    
    Coming in later releases: Files to process are specified for all (--process-all), 
    undone (--process_undone), or a custom list (--process_custom) of videos. If a custom
    list, include the path to a process file (.prc) containing all files to process, 
    each on a separate line. 
    
    Additional flags following '--batch' for concatenating ('--concat') all results 
    files and performing a repeated measures ANOVA statistics ('--stats') for longitudinal
    studies can also be specified.
    ----
    Inputs:
      None
    ----
    Returns:
      args (list): list of arguments passed to program
    '''
    parser = argparse.ArgumentParser(prog='FreeClimber',
                                    description='Calculates the climbing velocity of flies in a negative geotaxis assay',
                                    usage='%(prog)s [options] path',
                                    epilog='For documentation and a tutorial, see https://github.com/adamspierer/FreeClimber',
                                    allow_abbrev=False)

    ## Specifying file (video file for gui or configuration file for batch)
    parser.add_argument('--config_file', 
                        type=str, 
                        required=True, 
                        help="Path to configuration file (ends with '.cfg')")

    ## Insert 3
    ## Insert 8
    
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
    return args


## Checking argparse arguments
def check_args(args):
    '''Checks arguments passed to parser and kills program if invalid path.
    ----
    Inputs:
      args (list): Arguments passed to parsed
    ----
    Returns:
      file (str): File path if it is entered and valid, or None if not.
    '''

    if args.debug: print('-------->  check_args')
    ## Confirm file path is valid
    if os.path.isfile(args.config_file) and args.config_file.endswith('.cfg'):
        print('--> Configuraton file: ',args.config_file)
    elif args.config_file != None and ~os.path.isfile(args.config_file):
        print("--> Exiting program: Invalid file path, not a configuration file, or file lacks '.cfg' suffix")
        print('!Not valid!', args.config_file,'\n')
        raise SystemExit
    
    ## Insert 5
    ## Insert 2
    return args.config_file
    
    
def startup():
    '''Function prints top line and runs argument parsing function.
    ----
    Inputs:
      None
    ----
    Returns:
      args (list): list of arguments passed to program
    '''
    line_length = 72
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    line1 = '## FreeClimber v.%s ' % str(version)
    line2 = "## Beginning program @ %s " % str(now)

    print('\n'+'#' * line_length)
    print(line1 + '#'*(line_length-len(line1)))
    print(line2 + '#'*(line_length-len(line2)))
    print('#' * line_length)
    return         

def main():
    startup()
    args = define_argument_parser()
    if args.debug: print(args,'\n')
    cfg_file = check_args(args)
    variables = load_variables(config_file = cfg_file)
    global skip
    skip = False

    ## Set up FreeClimber object and load parameters
    fc = FreeClimber(cfg_file = cfg_file)
    fc.cfg_file = cfg_file
    fc.load_parameters()
    fc.setup()
    fc.file_list = fc.file_walker(fc.path_project, endswith=fc.file_suffix)

    ## Insert 1
    
    for File in fc.file_list:
        print(File)
        while skip == False:
            t0 = time()
            fc.count += 1
            fc.name = os.path.split(File)[-1]
            fc.process(File,variables)
            fc.timer(t0)
            skip = True
        skip = False

    ## insert 7
    fc.concat_slopes()
    fc.print_closing()
    return
        
if __name__ == '__main__':
    main()
    exit()
    
    
## Insert 1
#     fc.check_variables() # Not ready for this release

#     ## Process all or those unprocessed    
#     if fc.args.process_all or (~fc.args.process_all and ~fc.args.process_undone):
#         fc.file_list = fc.file_walker(fc.path_project, endswith=fc.file_suffix)
#     if fc.args.process_undone:
#         fc.file_list = fc.undone(fc.path_project)
#     fc.print_new_project()

## Insert 2
#     ## Checking which files to process
#     if args.process_all == None or args.process_all == True:
#         print('--> Processing all video files')
#     elif args.process_undone:
#         print('--> Processing undone video files')
#     elif args.process_custom != None:
#         print('--> Processing a list of files in a specified file a custom file list')
#         if os.path.isfile(args.process_custom) and args.process_custom.endswith('.prc'):
#             process_custom = self.read_custom(args.process_custom)
#         else:
#             print("--> Exiting program: Invalid file path or format, make sure proper file is specified.")
#             print('!!', args.process_custom)
#             raise SystemExit

## Insert 3
#     ## Batch processing: specifying which files to process
#     group_method = parser.add_mutually_exclusive_group(required=False)
#     group_method.add_argument('--process_all', 
#                         default=True, 
#                         action='store_true', 
#                         help="Process all files")
# 
#     group_method.add_argument('--process_undone', 
#                         default=False,
#                         action='store_true', 
#                         help="Process previously unprocessed files")
#     group_method.add_argument('--process_custom', 
#                         default=False, 
#                         type=str,
#                         help="Process custom list of files, must include file path as argument")

## Insert 4
#     ## Reading file with video paths for --process_custom argument
#     def read_custom(self, file):
#         '''Reads in process file (.prc) containing individual video file paths to process.
#         ----
#         Inputs:
#           file (str): Path to custom file (.prc)
#           suffix (str): Video file suffix read in from configuration (.cfg) file
#         ----
#         Returns:
#           _lines (list): List of files read from custom file (.prc)
#         '''
#         ## Read in lines from custom file
#         with open(file) as f:
#             lines = [line.rstrip() for line in f]
#         f.close()
#     
#         ## Check files are real
#         _lines = [vid for vid in lines if os.path.isfile(vid)]
#         print('----> %s of %s specified files have valid paths' % (len(lines),len(lines)))
#         return _lines

## Insert 5
#     ## Resolving conflicting commands for which to process (hierarchy: custom > undone > all)
#     if args.process_all:
#         if args.process_undone:
#             args.process_all = False
#         if args.process_custom != False:
#             args.process_all = False
#     if args.process_undone == True and args.process_custom != False:
#         args.process_undone = False

## Insert 6
#     def check_variables(self):
#         ## Making sure diameter is an odd number
#         if int(self.diameter) % 2 == False:
#             self.diameter = int(self.diameter) + 1
#             
#         ## Should probably check other variables at some point...
#         
#         return

## Insert 7
    # if args.concat:
#         fc.concat_slopes()
#     

## Insert 8
# ## Batch processing: additional arguments
#     parser.add_argument('--concat', 
#                         required=False, 
#                         default=True, 
#                         action='store_true',
#                         help="Set to False to prevent results files from concatenating")
