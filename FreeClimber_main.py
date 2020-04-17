#!/usr/bin/env python
# -*- coding: utf-8 -*- 

## Importing external packages
import os
import detector
import sys
from time import time
from numpy import sort, repeat
from matplotlib.cm import Greys_r
import pandas as pd
import matplotlib.pyplot as plt


def load_variables(config_file):
    ## Finding the path_project
    if os.path.isfile(config_file):
        pass
    else:
        print('\n\nExiting program. Invalid configuration file')
        raise SystemExit ## Issue when running GUI, terminates this script, but not the GUI

    with open(config_file,'r') as f:
        print('## Reading in configuration file:', config_file)
        variables = f.read().split('\n')
    f.close()

    ## Importing configuration file variables
    variables = [item for item in variables if ~item.startswith(' ') or item.startswith('#')]
    return variables

def check_arguments():
    ## I want to add arguments for which videos things to processes, so the GUI doesn't do all, or someone can specify 'how' they want it done: 1, not yet processed (undone), or all
    args = sys.argv[1:]
    if len(args) == 0 or len(args) > 2:
        print('\n\nExiting program. Too many or too few arguments -- required argument 1 is .cfg file path and optional argument 2 is a boolean to concatenate all slope files at the end')
        raise SystemExit ## Issue when running GUI, terminates this script, but not the GUI

    if os.path.isfile(args[0]) and args[0].endswith('.cfg'):
        config_file = args[0]
    else:
        print('\n\nExiting program. Invalid configuration file')
        raise SystemExit ## Issue when running GUI, terminates this script, but not the GUI
    try:
        if type(args[1]) == bool:
            pass
    except:
        args.append(False)
    return args


class free_climber(object):
    def __init__(self, cfg_file):
        self.cfg_file = cfg_file
        
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
        
    def file_walker(self,folder = None,endswith=None):
        '''
        From a specified parent folder, find all child files with the specified suffix
        '''
        _lst=[]
        for root, dirs, files in os.walk(folder):
            for name in files:
                if name.endswith(endswith):
                    _lst.append((os.path.join(root, name)))
        return sorted(_lst)

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
        self.file_list = self.file_walker(self.path_project,self.file_suffix)
        self.todo = self.undone(self.path_project)
        self.count = 0
        self.first_run = True
        self.print_new_project()
        return
        
    def print_new_video(self,name):
        print('='*72)
        line_to_print = '== [%s || %s] ' % (self.count,len(self.todo)), self.name
        line_to_print = ('').join(line_to_print)
        print(('').join(line_to_print), '=' * (71 - len(line_to_print)))
        print('='*72)
        return

    def print_new_project(self):
        print("##\n## Analyzing videos in: %s \n##" % self.path_project )
        print("## Of the %s videos in path_project, %s remain unprocessed\n##\n##" % (str(len(self.file_list)), str(len(self.todo))) )
        return

    def process(self,File,variables):
        self.print_new_video(self.name)

        d = detector.detector(File)
        d.step_1()#File)#,variables)
        d.step_2()    
#         d.step_3()
        d.step_4()
        ## -- [ Step 5 ] - Calculating local linear regressions
#         d.step_5()

        ## -- [ Step 6 ] - Plotting data
        d.step_6()

        ## -- [ Step 7 ] - Setting up slopes file
        d.step_7()
        
        self.first_run = False
        return

    def concat_slopes(self):
        print("\nFinal step: Concatenating slope files")
        slope_files = self.file_walker(self.path_project,endswith='.slopes.csv')
        to_concat = [slope_file for slope_file in slope_files]
        print("    - Concatenating",len(to_concat),"files")
        self.slopes = pd.concat([pd.read_csv(item) for item in to_concat])
        self.path_result = self.path_project+'/results.csv'
        self.slopes.to_csv(self.path_result,index=False)
        return
    
    def check_variables(self):
        ## Making sure diameter is an odd number
        if int(self.diameter) % 2 == False:
            self.diameter = int(self.diameter) + 1
            
        ## Should probably check other variables at some point...
        
        return

def print_opening():
    ## First few lines when running script
    print('\n')
    print('#' * 72)
    print('## Vial Detector 1.0','#'*51)
    print('#' * 72)
    return         

def main(cfg_file):
    print_opening()
    variables = load_variables(cfg_file)
    vd = free_climber(cfg_file = cfg_file)
    vd.cfg_file = cfg_file
    vd.load_parameters()
    vd.setup()
    vd.check_variables()
#     for File in vd.file_list:
    for File in vd.todo:
        t0 = time()
        vd.count += 1
        vd.name = os.path.split(File)[-1]
        vd.process(File,variables)
        vd.timer(t0)

    if args[1]:
        vd.concat_slopes()
    
    return
        
if __name__ == '__main__':
    args = check_arguments()
    cfg_file = args[0]
    main(cfg_file)
    exit()