#!/usr/bin/env python

## File name : gather_files.py
## Created by: Adam N. Spierer
## Date      : May 2020
## Purpose   : Used in conjunction with the FreeClimber software to search through a 
##              parent folder and extract all files with a common suffix

## Import modules
from numpy import unique
from os import walk,path
import argparse
from datetime import datetime

## FreeClimber version
version = '0.3.1'

def define_argument_parser():
    '''Defines arguments to be parsed, via argparse module.
    
    Takes in arguments for the file suffix (suffix) to look for in the parent folder 
      (parent_folder). Other arguments to search for files without a '.slopes.csv' 
      counterpart (undone) and two output methods: print to command line (print_files)
      and save the processed file (save_files; as custom.prc).
    ----
    Inputs:
      None
    ----
    Returns:
      args (object): Namespace object containing the flags and arguments passed to program
    '''    
    ## Set up the command line parser
    parser = argparse.ArgumentParser(prog='FreeClimber',
                                    description='Calculates the climbing velocity of flies in a negative geotaxis assay\ngather_files.py - Gathers all files with a common suffix in a parent folder and saves to a process (.prc) file',
                                    #  usage='%(prog)s [options] path',
                                    epilog='For documentation and a tutorial, see https://github.com/adamspierer/FreeClimber',
                                    allow_abbrev=False)

    ## Specifying configuration file
    parser.add_argument('--parent_folder', 
                        type=str, 
                        required=True, 
                        help="Path to parent ")
    
    parser.add_argument('--suffix',
                        type=str,
                        required=True,
                        help="Specify a common file suffix (h264, avi, mp3)")
    
    parser.add_argument('--undone', 
                        default=False,
                        action='store_true', 
                        help="Gather all undone files")

    parser.add_argument('--save_files', 
                        default=False,
                        action='store_true', 
                        help="Saves a processed (.prc) file to the parent_folder")

    parser.add_argument('--print_files', 
                        default=False,
                        action='store_true', 
                        help="Prints file list to command line")


    ## Final parser specifications
    parser.version = version ## Program version
    parser.add_argument('-v','--version', 
                        action='version')

    args = parser.parse_args()
    return args
    
def file_walker(folder = None, endswith = None, undone = None):
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
    _list1,_list2 = [],[]
    for root, dirs, files in walk(folder):
        for name in files:
            if name.endswith(endswith):
                _list1.append((path.join(root, name)))
            if undone:
                if name.endswith('.slopes.csv'):
                    _list2.append(path.join(root, name[:-11])+'.'+endswith)

    ## Return a sorted list of all files with the file suffix
    if undone == False:
        _list = sorted(unique(_list1))
        return _list

    ## Return a sorted list of all undone files            
    if undone:
        _list1,_list2 = set(_list1),set(_list2)
        _list = list(_list1.difference(_list2))
        _list = [item+'.'+endswith for item in _list]
        _list = sorted(_list)
        
        if ~len(_list):
            print('All files previously processed, re-evaluate your inputs if this message is a surprise.')
        return _list


def export(save_files = False, print_files = False, file_list=None, undone=False,
            suffix = None, destination = None):
    '''Exporting list of files to the command line, a new file, or both'''
    if print_files:
        ## Generating text if undone is true
        if undone: _str = 'previously unprocessed (no .slopes.csv file)'
        else: _str = ''
        
        ## Printing list header in command line output
        print("Exporting %s%s file(s) with suffix '%s' found in: %s" % (len(file_list), 
                                                                        _str,
                                                                        suffix,
                                                                        destination))
        print('------------')
        ## Printing each item found
        for item in file_list:
            print(item)
        print('\n')
    if save_files:
        
        save_path = destination + 'custom.prc'
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        print('Saving files to:', save_path)
        with open(save_path, 'w') as f:
            print('## FreeClimber ##', file=f)
            print('## Custom file list generated from gather_files.py ##', file=f)
            print('## Generated @ ' + now, file = f)
            print('##',file = f)

            for item in file_list:
                print(item,file = f)
        f.close()
    return

def main():
    '''Main function to run all sub-functions'''
    
    ## Define the command line interface arguments
    print('\n')
    args = define_argument_parser()
    
    ## Search through the parent folder, looking for files w/ suffix. Option for undone files
    if path.isdir(args.parent_folder): 
        file_list = file_walker(folder = args.parent_folder,
                                endswith = args.suffix,
                                undone = args.undone)
        print('\n')
    
    ## Exporting the identified files, to the command line and/or save as a file in parent
    export(save_files = args.save_files, 
           print_files = args.print_files, 
           file_list = file_list,
           suffix = args.suffix,
           destination = args.parent_folder)
    return
    
## Main function
if __name__ == '__main__':
    main()