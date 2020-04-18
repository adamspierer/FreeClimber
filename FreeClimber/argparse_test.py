#!/usr/bin/env python

import argparse
import os

def read_custom(file):
    '''Reads in process file (.prc) containing individual video file paths to process.
    ----
    Inputs:
      file (str): Path to custom file (.prc)
      suffix (str): Video file suffix read in from configuration (.cfg) file
    ----
    Returns:
      _lines (list): List of files read from custom file (.prc)
    '''
    ## Read in lines from custom file
    with open(file) as f:
        lines = [line.rstrip() for line in f]
    f.close()
    
    ## Check files are real
    _lines = [vid for vid in lines if os.path.isfile(vid)]
    print('----> %s of %s specified files have valid paths' % (len(lines),len(lines)))
    return _lines
    
    
def define_argument_parser():
    '''Defines arguments to be parsed, via argparse module. 
    Given a method argument, 'gui' or 'batch', different downstream arguments are required.
    
    For 'gui', an optional 'file' argument specifies which video to being with, though
    if not specified the program will open a browser to select the file. 
    
    For 'batch', the '--file' flag is required and points to a configuration file (.cfg)
    created from the GUI. Files to process are specified for all (--process-all), 
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

    ## Analysis platform
    group_method = parser.add_mutually_exclusive_group(required=True)
    group_method.add_argument('-b', '--batch', 
                                action='store_true', 
                                help="Run from command line")
    group_method.add_argument('-g', '--gui', 
                                action='store_true', 
                                help="Run from GUI")

    ## Specifying file (video file for gui or configuration file for batch)
    parser.add_argument('--file', 
                        type=str, 
                        required=False, 
                        help="Path to file (.cfg file if -b/--batch; (optional) video file if '-g/--gui')")

    ## Batch processing: specifying which files to process
    parser.add_argument('--process_all', 
                        default=True, 
                        action='store_true', 
                        help="With 'batch', process all files")
    parser.add_argument('--process_undone', 
                        default=False,
                        action='store_true', 
                        help="With 'batch', process previously unprocessed files")
    parser.add_argument('--process_custom', 
                        default=False, 
                        type=str,
                        help="With 'batch', process custom list of files. Include file path as argument")

    ## Batch processing: additional arguments
    parser.add_argument('--concat', 
                        required=False, 
                        default=True, 
                        action='store_true',
                        help="With 'batch', concatenate results into single file")
    parser.add_argument('--stats', 
                        required=False, 
                        default=False, 
                        action='store_true',
                        help="With 'batch', performs a repeated measures ANOVA. See documentation for caveats.")

    ## Final parser specifications
    parser.version = '1.0' ## Program version
    parser.add_argument('-v','--version', 
                        action='version')
    args = parser.parse_args()
    return args


## Check input arguments
def startup():
    '''Function prints top line and runs argument parsing function.
    ----
    Inputs:
      None
    ----
    Returns:
      args (list): list of arguments passed to program
    '''
    print('#'*80+'\n'+'Running FreeClimber\n')
    args = define_argument_parser()
    return args


def check_args(args):
    '''Checks arguments passed to parser.
    Inputs:
      args (list): Arguments passed to parsed
    ----
    Returns:
      None
    '''
    print('check_args')
    if args.gui:
        ## Confirm file path is valid
        if args.file != None and os.path.isfile(args.file):
            print('--> Video file: ',args.file)
        elif args.file != None and ~os.path.isfile(args.file):
            print('--> Exiting program: Invalid file path. Make sure video file for GUI is properly specified.')
            print('!!', args.file)
            raise SystemExit
        else:
            print('--> No video file specified, select from Browser')
            args.file = ''

    else:
        print('Batch')
        ## Checking configuration file
        if args.file == None:
            print('--> Exiting program: Must include path to configuration (.cfg) file')
            raise SystemExit
        elif os.path.isfile(args.file) and args.file.endswith('.cfg'):
            pass
        else:
            print('--> Exiting program: Invalid configuration file, make sure proper file with ".cfg" is specified.')
            print('!!', args.file)
            raise SystemExit
    
        ## Checking which files to process
        if args.process_all == None:
            print('--> Processing all video files')
        elif args.process_undone:
            print('--> Processing undone video files')
        elif args.process_custom != None:
            print('--> Processing a list of files in a specified file a custom file list')
            if os.path.isfile(args.process_custom) and args.process_custom.endswith('.prc'):
                process_custom = read_custom(args.process_custom)
            else:
                print("--> Exiting program: Invalid file path or format, make sure proper file is specified.")
                print('!!', args.process_custom)
                raise SystemExit
    return True


def main():
    '''Running main function'''
    args = startup()
    check_args(args)
    if args.gui:
        print('COMMAND: pythonw ./FreeClimber_gui.py %s' % args.file)
    if args.batch:
        print('COMMAND: python ./FreeClimber_batch.py %s' % args.file)
    print(args)
    print('Done!')
    return


if __name__ == "__main__":
    main()
