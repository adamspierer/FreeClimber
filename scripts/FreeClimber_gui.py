#!/usr/bin/env python
# -*- coding: utf-8 -*- 

## File name : FreeClimber_main.py
## Created by: Adam N. Spierer
## Date      : May 2020
## Purpose   : Graphical User Interface wrapper for FreeClimber

## Version number
version = '1.0.4'
doi =  'https://doi.org/10.1242/jeb.229377' ## Link to published paper

## More universal modules
import gc
import wx
import os
import sys
import time
import argparse
import traceback
import matplotlib
## Select the wx backend before pyplot is imported (anywhere, including by
## detector_fng) so pyplot does not bind to a different backend first.
matplotlib.use('WXAgg')
import numpy as np
import pandas as pd
from datetime import datetime

## Matplotlib configurations
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
import matplotlib.pyplot as plt

## Local imports
from detector_fng import detector


## Set up wxPython application window
class App(wx.App):
    '''Create the main window and insert the custom frame'''
    def OnInit(self):
        frame = create(None, file_name)
        self.SetTopWindow(frame)
        frame.Show(True)
        return True

## Main GUI frame
class main_gui(wx.Frame):    
    '''Setting up the placement and '''
    def __init__(self, parent, file_path):
        if args.debug: print('main_gui.__init__')
        self.args = args
        self.initialize_controls(parent)
        self.box_sizer.Add(self.panel1, 0, wx.EXPAND)

        ## Initialize GUI plot
        self.figure = Figure()
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.box_sizer.Add(self.canvas, 1, wx.EXPAND)

        ## Default rectangle variable
        self.pressed=False
        
        ## Initialize bottom text bar
        self.box_sizer.Add(self.status_bar, 0, border=0, flag=0)
        self.status_bar.SetStatusText('Ready', 0)
        rect = self.status_bar.GetFieldRect(0)
    
        ## Set up names
        self.video_file = file_path

        ## Create a dialog box at the beginning if the video path is not a real file
        if self.args.video_file == None:
            openFileDialog = wx.FileDialog(self, "Open Video file", "", "",
                                       "Video files (*.*)|*.*", 
                                       wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

            ## Closing program if Browse box is closed
            if openFileDialog.ShowModal() == wx.ID_CANCEL:
                print('\nExiting Program. Must select a video from box or enter path into command line\n')
                raise SystemExit

            ## Setting video name from dialog box selection
            self.video_file = openFileDialog.GetPath()        

        ## Passing individual inputs from GUI to a list
        self.update_names()
        self.input_names = ['x','y',
                            'w','h',
                            'check_frame','blank_0',
                            'blank_n','crop_0',
                            'crop_n','threshold',
                            'diameter','minmass',
                            'maxsize','ecc_low',
                            'ecc_high','vials',
                            'window','pixel_to_cm',
                            'frame_rate','vial_id_vars',
                            'outlier_TB','outlier_LR',
                            
                            'naming_convention','path_project','file_suffix',
                            
                            'convert_to_cm_sec','trim_outliers']
        self.parameter_names = ['self.input_'+item for item in self.input_names]
        self.input_values = [## int
                             self.input_x,self.input_y,
                             self.input_w,self.input_h,
                             self.input_check_frame,self.input_blank_0,
                             self.input_blank_n,self.input_crop_0,
                             self.input_crop_n,self.input_threshold,
                             self.input_diameter,self.input_minmass,
                             self.input_maxsize,self.input_ecc_low,
                             self.input_ecc_high,self.input_vials,
                             self.input_window,self.input_pixel_to_cm,
                             self.input_frame_rate,self.input_vial_id_vars,
                             self.input_outlier_TB,self.input_outlier_LR,
                             ## str
                             self.input_naming_convention,self.input_path_project,self.input_file_suffix,
                             ## bool
                             self.input_convert_to_cm_sec,self.input_checkBox_trim_outliers]

        ## Enable all buttons
        button_list = ['browse_video','reload_video','test_parameters','store_parameters']
        for button in button_list:
            getattr(self, 'button_' + button).Enable(True)

        ## Load video
        self.load_video()
        self.status_bar.SetStatusText("Ready...",0)
        if args.debug: print('End of init')
        return

    def update_variables(self):
        '''Updates the detection variables'''
        if args.debug: print('main_gui.update_variables')
        variables = []
        
        
        ## Including integers
        for item,jtem in zip(self.input_names[:22],self.input_values[:22]):
#             print('int',item,jtem)
            phrase = str(item + '='+ jtem.GetValue())
            if args.debug: print('    '+phrase)
            variables.append(phrase)

        ## Including strings - type I. Paths are written with forward slashes:
        ## the batch runner reads them back fine on Windows, and hand-edited
        ## copies of the .cfg stay free of backslash escapes.
        for item,jtem in zip(self.input_names[22:24],self.input_values[22:24]):
#             print('str-I',item,jtem)
            value = jtem.GetValue()
            if item == 'path_project': value = value.replace('\\', '/')
            phrase = str(item + '="'+ value +'"')
            if args.debug: print('    '+phrase)
            variables.append(phrase)

        ## Including strings - type II. file_suffix is read from the current
        ## video name each time (input_values holds a copy made at start-up,
        ## which went stale after browsing to a video with another extension).
        phrase = 'file_suffix="%s"' % self.input_file_suffix
        if args.debug: print(phrase)
        variables.append(phrase)
        
        ## Including booleans
        for item,jtem in zip(self.input_names[25:],self.input_values[25:]):
#             print('bool',item,jtem)
            phrase = str(item + '=%s' % str(jtem.GetValue()))
            if args.debug: print('    '+phrase)
            variables.append(phrase)

        ## Step 6: per-fly trajectory analysis (linking + tortuosity/meandering).
        ## Appended here (not via the index-sliced input_names/input_values lists)
        ## so the fixed slicing above is left untouched. Tortuosity is gated on
        ## individual mode in step_5, so checking it implies individual mode.
        tortuosity = self.input_tortuosity_enabled.GetValue()
        individual = self.input_analysis_individual.GetValue() or tortuosity
        variables.append("analysis_mode=%s" % ("'individual'" if individual else "'cohort'"))
        variables.append('tortuosity_enabled=%s' % str(tortuosity))
        ## Linking + tortuosity numeric parameters (documented defaults; fine-tune
        ## in the saved .cfg if needed). Written whenever the config is saved.
        variables.append('link_search_range=15')
        variables.append('link_memory=3')
        variables.append("link_predictor='nearest_velocity'")
        variables.append('link_min_track_length=5')
        variables.append('tortuosity_smoothing_window=5')
        variables.append('tortuosity_velocity_threshold=1.0')
        variables.append('tortuosity_bout_min_frames=10')
        variables.append('tortuosity_bout_min_displacement=2.0')
        if args.debug:
            print('    analysis_mode individual=%s, tortuosity_enabled=%s' % (individual, tortuosity))

        ## FNG detection parameters (documented defaults; fine-tune in the .cfg)
        variables.append('fng_enabled=True')
        variables.append('fng_smooth_window=5')
        variables.append('fng_climb_thresh=0.10')
        variables.append('fng_fall_thresh=0.10')
        variables.append('fng_min_gap=5')
        variables.append('fng_min_range_cm=2.0')

        ## Failure to climb (Step 7 row) and the vial floor
        variables.append('ftc_enabled=True')
        ## Blank line height -> None -> the line is the top of the drawn ROI box
        variables.append('ftc_height_cm=%s' % self._number_or_none(self.input_ftc_height_cm))
        variables.append('ftc_window_sec=%s' % self._number_or_none(self.input_ftc_window_sec))
        variables.append('ftc_eval_frames=5')
        variables.append('ftc_min_coverage=0.8')
        variables.append('flies_per_vial=%s' % self._flies_per_vial())
        floor = self._number_or_none(self.input_floor_y)
        if floor is not None:
            ## Entered in full-image pixels (as read off the video axes);
            ## the detector wants pixels from the top of the ROI.
            try: floor = floor - int(self.input_y.GetValue())
            except ValueError: pass
        variables.append('floor_y=%s' % floor)
        return variables

    @staticmethod
    def _number_or_none(control, default=None):
        '''Numeric value of a text box, or default when blank/invalid.'''
        text = control.GetValue().strip()
        if text == '':
            return default
        try:
            value = float(text)
        except ValueError:
            return default
        return int(value) if value.is_integer() else value

    def _flies_per_vial(self):
        '''Flies loaded per vial: blank -> None, '10' -> 10, '10,10,9' -> [10, 10, 9].'''
        tokens = [t.strip() for t in self.input_flies_per_vial.GetValue().split(',') if t.strip()]
        try:
            counts = [int(t) for t in tokens]
        except ValueError:
            return None
        if not counts:
            return None
        return counts[0] if len(counts) == 1 else counts

    def load_video(self):
        '''Function for loading the video when the respective button is pressed'''
        if args.debug: print('main_gui.load_video')
        
        ## Set up
        self.status_bar.SetStatusText("Loading video",0)
        self.figure.clear()
        self.axes   = [self.figure.add_subplot(111), ]
        
        ## Confirm file is a path, or folder has specified suffix
        status = self.check_specified_video()
        if status:
            print("Loading:", self.video_file)
            self.update_names()
            for item in self.input_values[:4]:
                item.SetEditable(True)
                item.Enable(True)
            self.checkBox_fixed_ROI.Enable(True)
            self.input_convert_to_cm_sec.Enable(True)

            ## Busy cursor while the detector object is called and initialized.
            ## A video that cannot be read raises SystemExit in the detector;
            ## report it instead of letting it close the GUI.
            wx.BeginBusyCursor()
            try:
                vars = self.update_variables()
                self.detector = detector(self.video_file,
                                        gui=True,
                                        variables = vars)

                self.axes[0].imshow(self.detector.image_stack[0])
                self.figure.canvas.draw()
            except (Exception, SystemExit) as e:
                self.detector = None
                self.report_error('Could not load the video', e)
                return
            finally:
                wx.EndBusyCursor()

            ## Setting mechanism for drawing the ROI rectangle. Connect the
            ## handlers once: reloading used to stack duplicate callbacks.
            if not getattr(self, 'handlers_connected', False):
                self.canvas.Bind(wx.EVT_ENTER_WINDOW, self.ChangeCursor)
                self.canvas.mpl_connect('button_press_event', self.draw_rectangle)
                self.canvas.mpl_connect('button_release_event', self.on_release)
                self.canvas.mpl_connect('motion_notify_event', self.on_motion)
                self.handlers_connected = True
            self.rect = Rectangle((0,0), 1, 1, fill=False, ec='r')
            self.axes[0].add_patch(self.rect)

            ## Auto-set GUI parameters from the video
            n_frames = self.detector.n_frames
            self.input_blank_0.SetValue('0')
            self.input_blank_n.SetValue(str(n_frames))
            self.input_crop_0.SetValue('0')
            self.input_crop_n.SetValue(str(n_frames))
            self.input_ecc_low.SetValue('0')
            self.input_ecc_high.SetValue('1')
            self.input_path_project.SetValue(self.folder.replace('\\', '/'))
            self.input_naming_convention.SetValue(self.name)
            self.input_vial_id_vars.SetValue(str(len(self.input_naming_convention.GetValue().split('_'))))

            ## Frame rate comes from the video itself (a wrong frame rate skews
            ## every velocity); fall back to the text box if it is unavailable.
            fps = getattr(self.detector, 'video_fps', None)
            if fps:
                self.input_frame_rate.SetValue(str(int(round(fps))) if abs(fps - round(fps)) < 0.01
                                               else '%.3f' % fps)
            try:
                fps = float(self.input_frame_rate.GetValue())
            except ValueError:
                fps = 30.0

            ## Check frame: about t = 2 seconds, or the last frame of a short video
            self.input_check_frame.SetValue(str(int(min(fps * 2, n_frames - 1))))

            ## Try to make the local linear regression window size 2 seconds, but if not then 35% of the frames in the video
            if n_frames < fps*2:
                self.input_window.SetValue(str(int(n_frames * .35)))
            else:
                self.input_window.SetValue(str(int(fps*2)))

            ## Enable Test parameter button if disabled from prior testing
            self.button_test_parameters.Enable(True)
            self.x0,self.y0 = 0, 0
            self.x1,self.y1 = self.detector.width,self.detector.height
        
            ## Display the first frame of the video in the GUI
            self.update_ROIdisp()
            self.canvas.draw()
            self.status_bar.SetStatusText(
                'Loaded %s frames at %s fps. Draw the ROI with its bottom edge at or just '
                'below the vial floor.' % (n_frames, self.input_frame_rate.GetValue()), 0)
        else:
            return

    def report_error(self, what, error):
        '''Show an error in the console, the status bar and a message box.'''
        if isinstance(error, SystemExit):
            detail = 'See the console for details.'
        else:
            detail = '%s: %s' % (type(error).__name__, error)
            traceback.print_exc()
        print('!! %s. %s' % (what, detail))
        self.status_bar.SetStatusText('%s. %s' % (what, detail), 0)
        wx.MessageBox('%s.\n\n%s' % (what, detail), 'FreeClimber', wx.OK | wx.ICON_ERROR, self)

    def update_names(self):
        '''Updates the names of variables within the program. Generally variables set for naming files.'''
        if args.debug: print('main_gui.update_names')
        self.status_bar.SetStatusText("Updating file names...",0)
        self.text_video_path.SetLabelText(self.video_file)
        ## splitext copes with names containing extra dots (or none)
        self.folder, name = os.path.split(self.video_file)
        self.name, suffix = os.path.splitext(name)
        self.input_file_suffix = suffix.lstrip('.')

        ## Naming files to be generated
        self.name_noext  = os.path.join(self.folder,self.name)
        if args.debug: print('name:',self.name_noext,"+ file suffixes")

        ## Set path_project default to the folder of the selected video file
        if self.input_path_project.GetValue() == '':
            self.input_path_project.SetValue(self.folder.replace('\\', '/'))
        return
        
    def check_specified_video(self):
        if args.debug: print('main_gui.check_specified_video')
        self.status_bar.SetStatusText("Checking specified video...",0)
    
        ## Check file path and update names
        if os.path.isfile(self.video_file):
            self.button_reload_video.Enable(True)
            self.button_test_parameters.Enable(True)
            self.update_names()
            return True
        
        else:
            self.video_file = "No or invalid file entered. Please change the file path"
            self.button_browse_video.Enable(True)
            return False

    ## Commands for drawing the ROI rectangle
    def ChangeCursor(self, event):
        '''Change cursor into crosshair type when enter the plot area'''
        self.canvas.SetCursor(wx.Cursor(wx.CURSOR_CROSS))
        return

    def draw_rectangle(self, event):
        '''Draw ROI rectangle'''
        self.status_bar.SetStatusText("Drag to draw the ROI rectangle",0)
        self.pressed = True
        if self.checkBox_fixed_ROI.IsEnabled() and event.xdata is not None:
            try:
                self.x0 = int(event.xdata)
                self.y0 = int(event.ydata)
                
                ## If the fixed_ROI box is checked, handle values differently
                if self.checkBox_fixed_ROI.GetValue():
                    try:
                        self.x1 = self.x0 + int(self.input_w.GetValue())
                        self.y1 = self.y0 + int(self.input_h.GetValue())
                    except ValueError:
                        pass
                    self.rect.set_width(self.x1 - self.x0)
                    self.rect.set_height(self.y1 - self.y0)
                    self.rect.set_xy((self.x0, self.y0))
                    self.canvas.draw()
                
                ## Set the values in the GUI and program to drawn rectangle
                self.update_ROIdisp()
            except (TypeError, ValueError):
                pass
        return

    def on_release(self, event):
        '''When mouse is on plot and button is released, redraw ROI rectangle, update ROI values'''
        self.status_bar.SetStatusText("Specify the detector parameters...",0)
        self.pressed = False
        if self.checkBox_fixed_ROI.IsEnabled():
            if self.checkBox_fixed_ROI.GetValue():
                pass
            else:
                self.redraw_rect(event)
            self.update_ROIdisp()
        return

    def on_motion(self, event):
        '''If the mouse is on plot and if the mouse button is pressed, redraw ROI rectangle'''
        if self.pressed and self.checkBox_fixed_ROI.IsEnabled() and (not self.checkBox_fixed_ROI.GetValue()):
            # Redraw the rectangle
            self.redraw_rect(event)
            self.update_ROIdisp()
        return

    def redraw_rect(self, event):
        '''Draw the ROI rectangle overlay'''
        try:
            x1 = int(event.xdata)
            y1 = int(event.ydata)
            if any([self.x1!=x1, self.y1!=y1]):                
                self.x1 = x1
                self.y1 = y1
                self.rect.set_xy((self.x0, self.y0))
                self.rect.set_width(self.x1 - self.x0)
                self.rect.set_height(self.y1 - self.y0)
                
                self.canvas.draw()
            else:
                pass
        except (TypeError, ValueError):  # pointer outside the image (xdata None)
            pass
        return

    def update_ROIdisp(self):
        '''Updates the ROI coordinates as the rectangle is drawn. The rectangle
        may be dragged in any direction and past the image edge, so the values
        are normalised to the top-left corner and clipped to the image.'''
        width = getattr(getattr(self, 'detector', None), 'width', None)
        height = getattr(getattr(self, 'detector', None), 'height', None)
        x0, x1 = sorted((int(self.x0), int(self.x1)))
        y0, y1 = sorted((int(self.y0), int(self.y1)))
        x0, y0 = max(0, x0), max(0, y0)
        if width: x1 = min(x1, width)
        if height: y1 = min(y1, height)
        self.input_x.SetValue(str(x0))
        self.input_y.SetValue(str(y0))
        self.input_h.SetValue(str(max(0, y1 - y0)))
        self.input_w.SetValue(str(max(0, x1 - x0)))
        return

    def OnButton_testParButton(self, event):
        '''Tests the entered parameters when the `Test parameters` button is pressed'''
        if args.debug: print('main_gui.OnButton_testParButton')
        self.status_bar.SetStatusText("Testing parameters...",0)

        if getattr(self, 'detector', None) is None:
            self.report_error('No video loaded', ValueError('Browse to a video first'))
            return

        #Prep the parameters
        variables = self.update_variables()
        self.checkBox_fixed_ROI.Enable(False)

        ## Set up figure for plots — clear() removes old axes/artists; gc.collect()
        ## encourages immediate release of image arrays those artists held.
        self.figure.clear()
        gc.collect()
        self.axes = [self.figure.add_subplot(231),
                     self.figure.add_subplot(232),
                     self.figure.add_subplot(233),
                     self.figure.add_subplot(234),
                     self.figure.add_subplot(235),
                     self.figure.add_subplot(236)]

        ## Busy cursor while the main function runs. Detection problems (e.g.
        ## 'no spots post-filtering') raise SystemExit in the detector; show
        ## them instead of silently leaving blank plots.
        wx.BeginBusyCursor()
        try:
            variables = variables + ['debug='+str(args.debug)]
            self.detector.parameter_testing(variables, self.axes)
        except (Exception, SystemExit) as e:
            wx.EndBusyCursor()
            self.checkBox_fixed_ROI.Enable(True)
            self.figure.canvas.draw()
            self.report_error('Parameter test failed', e)
            return
        wx.EndBusyCursor()

        ## Renders plots in the GUI
        self.figure.tight_layout()
        self.figure.canvas.draw()
        
        # Enable buttons and print statements once parameter testing is complete
        self.button_reload_video.Enable(True)
        self.button_store_parameters.Enable(True)
        if args.debug: print('Parameter testing complete')
        self.status_bar.SetStatusText("Refine detector parameters by reloading the video, or finish optimization by pressing 'Save configuration'",0)
        return

    def OnButton_strParButton(self, event):
        '''Runs the 'save_parameter' function for creating the configuration file'''
        if args.debug: print('main_gui.OnButton_strParButton')
        self.save_parameter()
        self.button_store_parameters.SetBackgroundColour(wx.Colour(241,241,241))

    def set_config_file(self):
        '''Path of the configuration file: <project folder>/<video name>.cfg, or
        beside the video when the project path is not an existing folder.'''
        if args.debug: print('main_gui.set_config_file')
        ## Read the text box each time (it used to be replaced by a string on
        ## the first save, so later edits to the project path were ignored)
        project = self.input_path_project.GetValue().strip()
        if os.path.isdir(project):
            self.path_parameters = os.path.join(project, self.name + '.cfg')
        else:
            self.path_parameters = self.name_noext + '.cfg'
        return self.path_parameters

    def save_parameter(self):
        '''
        Save parameters as python list.
        New parameter sets appended to the config file.
        Each parameter sets come with a comment line, contain the datetime of analysis
        '''
        if args.debug: print('main_gui.save_parameter')

        variables = self.update_variables()
        self.path_parameters = self.set_config_file()

        ## Printing output to configuration file
        print('Saving parameters to:', self.path_parameters)
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        try:
            with open(self.path_parameters, 'w') as f:
                print('## FreeClimber ##', file=f)
                print('## Generated from file: '+self.video_file,file = f)
                print('##     @ ' + now, file = f)
                print('##',file = f)
                print('## Analysis parameters:',file = f)
                for item in variables:
                    print(item,file = f)
        except OSError as e:
            self.report_error('Could not save the configuration', e)
            return
        print("Configuration settings saved")
        self.status_bar.SetStatusText('Configuration saved to %s' % self.path_parameters, 0)
        return

    def OnButton_Browse(self, event):
        if args.debug: print('main_gui.OnButton_Browse')
        openFileDialog = wx.FileDialog(self, "Open Video file", "", "",
                                       "Video files (*.*)|*.*", 
                                       wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        if openFileDialog.ShowModal() == wx.ID_CANCEL:
            pass
        else:
            ## load_video draws the first frame; the figure used to be cleared
            ## afterwards, leaving a blank canvas after every Browse.
            self.video_file = openFileDialog.GetPath()
            self.update_names()
            self.load_video()
        return

    def OnButton_LoadVideo(self,event):
        '''Calls function to load the video when the `reload` button is pressed'''
        if args.debug: print('main_gui.OnButton_LoadVideo')
        self.load_video()
        return
        
    def initialize_sizers(self):
        '''Initializes the GUI window orientation'''
        if args.debug: print('main_gui.initialize_sizers')
        # Generated method, do not edit
        self.box_sizer = wx.BoxSizer(orient=wx.VERTICAL)
        self.SetSizer(self.box_sizer)
        return

    def initialize_controls(self, prnt):
        '''Initializes the GUI controls, for which there are many'''
        if args.debug: print('main_gui.initialize_controls')
        # Generated method, do not edit
        wx.Frame.__init__(self, id=wxID_text_title, name='', parent=prnt,
              pos=wx.Point(100, 30), size=wx.Size(950, 795),
              style=wx.DEFAULT_FRAME_STYLE,
              title='FreeClimber')
        self.SetClientSize(wx.Size(950, 773))

        ######
        ## Inputs for ROI Rectangle
        self.panel1 = wx.Panel(id=wxID_panel_1, name='panel1',
              parent=self, pos=wx.Point(0, 0), size=wx.Size(950, 266),
              style=wx.TAB_TRAVERSAL)
       
        ## Step 1 boxes
        self.text_step_1a = wx.StaticText(id=wxID_text_step_1a,
              label=u'Step 1a: Specify a video', name='text_step_1a', parent=self.panel1, 
              pos=wx.Point(col1,10),
              size=wx.Size(box_dimensions),style = wx.ALIGN_CENTER)
                      
        ## Browse
        self.button_browse_video = wx.Button(id=wxID_browse_video,
              label=u'Browse...', name=u'button_browse_video', parent=self.panel1,
              pos=wx.Point(col1, 30), size=wx.Size(box_dimensions), style=0)
        self.button_browse_video.Bind(wx.EVT_BUTTON, self.OnButton_Browse,
              id=wxID_browse_video)
        
        self.text_step_1b = wx.StaticText(id=wxID_text_step_1b,
              label=u'Step 1b: Define options', name='text_step_1b', parent=self.panel1, 
              pos=wx.Point(col1,65),
              size=wx.Size(box_dimensions),style = wx.ALIGN_CENTER)

        ## Pixel to cm
        self.text_pixel_to_cm = wx.StaticText(id=wxID_text_pixel_to_cm,
              label=u"Pixels / cm:", name='text_pixel_to_cm',
              parent=self.panel1, pos=wx.Point(col1, 85), size=wx.Size(medium_box_dimensions),
              style=0)
        self.input_pixel_to_cm = wx.TextCtrl(id=wxID_input_pixel_to_cm,
              name=u'input_pixel_to_cm', parent=self.panel1, pos=wx.Point(col1 + 95, 85), 
              size=wx.Size(medium_box_dimensions), style=0, value=u"1")
        
        ## Frame Rate    
        self.text_frame_rate = wx.StaticText(id=wxID_frame_rate,
              label=u'Frames / sec:', name='text_frame_rate',
              parent=self.panel1, pos=wx.Point(col1,115), size=wx.Size(box_dimensions),
              style=0)
        self.input_frame_rate = wx.TextCtrl(id=wxID_frame_rate,
              name=u'input_frame_rate', parent=self.panel1, pos=wx.Point(col1 + 95,115),
               size=wx.Size(medium_box_dimensions), style=0, value='25')
        
        ## Check box to convert final slope to cm
        self.input_convert_to_cm_sec = wx.CheckBox(id=wxID_input_convert_to_cm_sec,
              label=u'Convert to cm / sec', name=u'input_convert_to_cm_sec',
              parent=self.panel1, pos=wx.Point(col1, 145), size=wx.Size(250,22),
              style=0)

        ## Step 2 boxes
        self.text_step_2 = wx.StaticText(id=wxID_text_step_2,
              label=u'Step 2: Select ROI', name='text_step_2', parent=self.panel1, 
              pos=wx.Point(col2,10),
              size=wx.Size(box_dimensions),style = wx.ALIGN_LEFT)
        ## X
        self.text_x = wx.StaticText(id=wxID_text_x,
              label=u'x-pos.', name='text_x', parent=self.panel1, 
              pos=wx.Point(col2, 30),
              size=wx.Size(medium_box_dimensions),style = wx.ALIGN_LEFT)
        self.input_x = wx.TextCtrl(id=wxID_input_x,
              name=u'input_x', parent=self.panel1, 
              style=0, value=u'0',
              pos=wx.Point(col2 + 55, 30),size=wx.Size(medium_box_dimensions))
        
        ## Y      
        self.text_y = wx.StaticText(id=wxID_text_y,
              label=u'y-pos.', name='text_y', parent=self.panel1,
              pos=wx.Point(col2, 55), size=wx.Size(medium_box_dimensions), style=wx.ALIGN_LEFT)
        self.input_y = wx.TextCtrl(id=wxID_input_y,
              name=u'input_y', parent=self.panel1, pos=wx.Point(col2 + 55,55),
              size=wx.Size(medium_box_dimensions), style=0, value=u'0')
        
        ## Width
        self.text_w = wx.StaticText(id=wxID_text_w,
              label=u'Width:', name='text_w', parent=self.panel1,
              pos=wx.Point(col2, 80), size=wx.Size(medium_box_dimensions), style=wx.ALIGN_LEFT)
        self.input_w = wx.TextCtrl(id=wxID_input_w,
              name=u'input_w', parent=self.panel1, pos=wx.Point(col2 + 55, 80),
              size=wx.Size(medium_box_dimensions), style=0, value=u'0')
        
        ## Height
        self.text_h = wx.StaticText(id=wxID_text_h,
              label=u'Height:', name='text_h', parent=self.panel1,
              pos=wx.Point(col2,105), size=wx.Size(medium_box_dimensions), style=wx.ALIGN_LEFT)
        self.input_h = wx.TextCtrl(id=wxID_input_h,
              name=u'input_h', parent=self.panel1, pos=wx.Point(col2 + 55, 105),
               size=wx.Size(medium_box_dimensions), style=0, value=u'0')
        
        ## ROI rectangle stays same dimensions but can be redrawn. Not critical to keep
        self.checkBox_fixed_ROI = wx.CheckBox(id=wxID_check_box_ROI,
              label=u'Fixed ROI Size?', name=u'checkBox_fixed_ROI',
              parent=self.panel1, pos=wx.Point(col2, 145), size=wx.Size(250,22),
              style=0)
        self.checkBox_fixed_ROI.SetValue(False)
        ######
        
       ## Detection parameters
        self.text_step_3 = wx.StaticText(id=wxID_text_step_3,
              label=u'Step 3: Specify spot parameters', name='text_step_3', parent=self.panel1, 
              pos=wx.Point(col3,10),
              size=wx.Size(100,22),style = wx.ALIGN_LEFT)
                      
        ## Expected spot diameter
        self.text_diameter = wx.StaticText(id=wxID_text_diameter,
              label=u'Diameter:', name='text_diameter', parent=self.panel1,
              pos=wx.Point(col3 ,30), size=wx.Size(medium_box_dimensions), style=0)
        self.input_diameter = wx.TextCtrl(id=wxID_input_diameter,
              name=u'input_diameter', parent=self.panel1, pos=wx.Point(col3 + 100,30), 
              size=wx.Size(medium_box_dimensions), style=0, value=u'7')
        
        ## Maximum spot diameter
        self.text_maxsize = wx.StaticText(id=wxID_text_maxsize,
              label=u'MaxDiameter:', name='text_maxsize', parent=self.panel1,
              pos=wx.Point(col3,55), size=wx.Size(medium_box_dimensions), style=0)
        self.input_maxsize = wx.TextCtrl(id=wxID_input_maxsize,
              name=u'input_maxsize', parent=self.panel1, pos=wx.Point(col3 + 100,55), 
              size=wx.Size(medium_box_dimensions), style=0, value=u'11')

        ## Minimum spot 'mass'
        self.text_minmass = wx.StaticText(id=wxID_text_minmass,
              label=u'MinMass:', name='text_minmass', parent=self.panel1,
              pos=wx.Point(col3,80), size=wx.Size(medium_box_dimensions), style=0)
        self.input_minmass = wx.TextCtrl(id=wxID_input_minmass,
              name=u'input_minmass', parent=self.panel1, pos=wx.Point(col3 + 100,80), 
              size=wx.Size(medium_box_dimensions), style=0, value=u'100')

        ## Spot threshold
        self.text_threshold = wx.StaticText(id=wxID_text_threshold,
              label=u'Threshold:', name='text_threshold', parent=self.panel1,
              pos=wx.Point(col3, 105), size=wx.Size(medium_box_dimensions), style=0)
        self.input_threshold = wx.TextCtrl(id=wxID_input_threshold,
              name=u'input_threshold', parent=self.panel1, pos=wx.Point(col3 + 100,105),
              size=wx.Size(medium_box_dimensions), style=0, value=u'"auto"')
        
        ## Eccentricity range
        self.text_ecc = wx.StaticText(id=wxID_text_ecc,
              label=u'Ecc/circularity:', name='text_ecc',
              parent=self.panel1, pos=wx.Point(col3, 130), size=wx.Size(medium_box_dimensions),
              style=0)
        self.input_ecc_low = wx.TextCtrl(id=wxID_input_ecc_low,
              name=u'input_ecc_low', parent=self.panel1, pos=wx.Point(col3+ 100, 130),
            size=wx.Size(small_box_dimensions), style=0, value=u'0')
        self.input_ecc_high = wx.TextCtrl(id=wxID_input_ecc_high,
              name=u'input_ecc_high', parent=self.panel1, pos=wx.Point(col3 + 140,130),
              size=wx.Size(small_box_dimensions), style=0, value=u'0')        
        
        #### Step 4 arguments
        ## Check frames
        self.text_step_4 = wx.StaticText(id=wxID_text_step_4,
              label=u'Step 4: Additional parameters', name='text_step_4', parent=self.panel1, 
              pos=wx.Point(col4,10),
              size=wx.Size(100,22),style = wx.ALIGN_LEFT)
        
        ## Background frames
        self.text_background_frames = wx.StaticText(id=wxID_text_background_frames,
              label=u'Background frames:', name='text_background_frames',
              parent=self.panel1, pos=wx.Point(col4, 30), size=wx.Size(medium_box_dimensions),
              style=0)
        self.input_blank_0 = wx.TextCtrl(id=wxID_input_blank_0,
              name=u'input_blank_0', parent=self.panel1, pos=wx.Point(col4+ 130, 30),
            size=wx.Size(small_box_dimensions), style=0, value=u'0')
        self.input_blank_n = wx.TextCtrl(id=wxID_input_blank_n,
              name=u'input_blank_n', parent=self.panel1, pos=wx.Point(col4 + 170,30),
              size=wx.Size(small_box_dimensions), style=0, value=u'0')

        ## crop frames
        self.text_crop_frames = wx.StaticText(id=wxID_text_crop_frames,
              label=u'Crop frames:', name='text_crop_frames',
              parent=self.panel1, pos=wx.Point(col4, 55), size=wx.Size(medium_box_dimensions),
              style=0)
        self.input_crop_0 = wx.TextCtrl(id=wxID_input_crop_0,
              name=u'input_crop_0', parent=self.panel1, pos=wx.Point(col4+ 130, 55),
            size=wx.Size(small_box_dimensions), style=0, value=u'0')
        self.input_crop_n = wx.TextCtrl(id=wxID_input_crop_n,
              name=u'input_crop_n', parent=self.panel1, pos=wx.Point(col4 + 170,55),
              size=wx.Size(small_box_dimensions), style=0, value=u'0')

        ## Check frames
        self.text_check_frames = wx.StaticText(id=wxID_text_check_frames,
              label=u'Check frame:', name='text_check_frames', parent=self.panel1,
              pos=wx.Point(col4, 80), size=wx.Size(115, 17), style=0)

        self.input_check_frame = wx.TextCtrl(id=wxID_input_check_frame,
              name=u'input_check_frame', parent=self.panel1, pos=wx.Point(col4 + 130, 80),
               size=wx.Size(small_box_dimensions), style=0, value=u'0')
        
        ## Vials
        self.text_vials = wx.StaticText(id=wxID_text_vials,
              label=u'Number of vials:', name='text_vials',
              parent=self.panel1, pos=wx.Point(col4, 105), size=wx.Size(133,22),
              style=0)
        self.input_vials = wx.TextCtrl(id=wxID_input_vials,
              name=u'input_vials', parent=self.panel1, pos=wx.Point(col4+130,105), 
              size=wx.Size(small_box_dimensions), style=0, value=u'1')
        
        ## Window size
        self.text_window = wx.StaticText(id=wxID_text_window,
              label=u'Window size:', name='text_window',
              parent=self.panel1, pos=wx.Point(col4, 130), size=wx.Size(133,22),
              style=0)
        self.input_window = wx.TextCtrl(id=wxID_input_window,
              name=u'input_window', parent=self.panel1, pos=wx.Point(col4+130,130), 
              size=wx.Size(small_box_dimensions), style=0, value='1')    
          

        ## Edge trim
        self.input_checkBox_trim_outliers = wx.CheckBox(id=wxID_check_box_outlier,
              label=u'Trim outliers? (TB                      LR)', name=u'checkBox_outlier',
              parent=self.panel1, pos=wx.Point(col4, 155), size=wx.Size(250,22),
              style=0)
        self.input_checkBox_trim_outliers.SetValue(False)
        self.input_outlier_TB = wx.TextCtrl(id=wxID_outlier_TB,
              name=u'input_outlier_TB', parent=self.panel1, pos=wx.Point(col4+130,155),
            size=wx.Size(small_box_dimensions), style=0, value=u'1')
        self.input_outlier_LR = wx.TextCtrl(id=wxID_outlier_LR,
              name=u'input_outlier_LR', parent=self.panel1, pos=wx.Point(col4+170,155),
              size=wx.Size(small_box_dimensions), style=0, value=u'3')


        self.text_step_5 = wx.StaticText(id=wxID_text_step_5,
              label=u'Step 5: Naming parameters', name='text_step_5', parent=self.panel1, 
              pos=wx.Point(col5,10),
              size=wx.Size(100,22),style = wx.ALIGN_LEFT)

        ## Naming convention
        self.text_naming_convention = wx.StaticText(id=wxID_text_naming_convention,
              label=u"Naming pattern:", name='text_naming_convention',
              parent=self.panel1, pos=wx.Point(col5, 30), size=wx.Size(medium_box_dimensions),
              style=0)
        self.input_naming_convention = wx.TextCtrl(id=wxID_input_naming_convention,
              name=u'input_naming_convention', parent=self.panel1, pos=wx.Point(col5,50), 
              size=wx.Size(large_box_dimensions), style=0, value='') 
        
        ## Variables
        self.text_vial_id_vars = wx.StaticText(id=wxID_text_vial_id_vars,
              label=u'Vial_ID variables:', name='text_vial_id_vars',
              parent=self.panel1, pos=wx.Point(col5, 80), size=wx.Size(133, 22),
              style=0)
        self.input_vial_id_vars = wx.TextCtrl(id=wxID_input_vial_id_vars,
              name=u'input_vial_id_vars', parent=self.panel1, pos=wx.Point(col5 + 130, 80),
              size=wx.Size(small_box_dimensions), style=0, value=u'2')

        ## Step 6: per-fly trajectory analysis (linking + tortuosity/meandering).
        ## These two checkboxes write the analysis_mode / tortuosity_enabled keys
        ## into the saved .cfg so the detector runs link_trajectories() and
        ## compute_tortuosity() (both gated on individual mode in step_5).
        self.text_step_6 = wx.StaticText(id=wx.ID_ANY,
              label=u'Step 6: Trajectory analysis', name='text_step_6', parent=self.panel1,
              pos=wx.Point(col5, 110),
              size=wx.Size(190, 22), style=wx.ALIGN_LEFT)
        self.input_analysis_individual = wx.CheckBox(id=wx.ID_ANY,
              label=u'Individual mode (per-fly linking)', name=u'input_analysis_individual',
              parent=self.panel1, pos=wx.Point(col5, 132), size=wx.Size(195, 22),
              style=0)
        self.input_analysis_individual.SetValue(False)
        self.input_tortuosity_enabled = wx.CheckBox(id=wx.ID_ANY,
              label=u'Tortuosity / meandering metrics', name=u'input_tortuosity_enabled',
              parent=self.panel1, pos=wx.Point(col5, 154), size=wx.Size(195, 22),
              style=0)
        self.input_tortuosity_enabled.SetValue(False)

        self.text_path_project = wx.StaticText(id=wxID_text_path_project,
              label=u"Project path:", name='text_path_project',
              parent=self.panel1, pos=wx.Point(col4-45, 180), size=wx.Size(medium_box_dimensions),
              style=0)
        self.input_path_project = wx.TextCtrl(id=wxID_input_path_project,
              name=u'input_path_project', parent=self.panel1, pos=wx.Point(col4 + 40, 180), 
              size=wx.Size(350,22), style=0, value='') 


        ## Step 7: failure to climb (and the vial floor it is measured from)
        row7 = 210
        self.text_step_7 = wx.StaticText(id=wx.ID_ANY,
              label=u'Step 7: Failure to climb', name='text_step_7', parent=self.panel1,
              pos=wx.Point(col1, row7 + 3), size=wx.Size(150, 22), style=0)
        self.text_ftc_height_cm = wx.StaticText(id=wx.ID_ANY,
              label=u'Line (cm):', name='text_ftc_height_cm', parent=self.panel1,
              pos=wx.Point(165, row7 + 3), size=wx.Size(60, 22), style=0)
        self.input_ftc_height_cm = wx.TextCtrl(id=wx.ID_ANY,
              name=u'input_ftc_height_cm', parent=self.panel1, pos=wx.Point(228, row7),
              size=wx.Size(small_box_dimensions), style=0, value=u'')
        self.input_ftc_height_cm.SetToolTip(
              'A fly that never gets this high above the vial floor within the time '
              'limit (and never falls) failed to climb. Blank = the top of the drawn '
              'ROI box.')
        self.text_ftc_window_sec = wx.StaticText(id=wx.ID_ANY,
              label=u'Within (s):', name='text_ftc_window_sec', parent=self.panel1,
              pos=wx.Point(275, row7 + 3), size=wx.Size(62, 22), style=0)
        self.input_ftc_window_sec = wx.TextCtrl(id=wx.ID_ANY,
              name=u'input_ftc_window_sec', parent=self.panel1, pos=wx.Point(340, row7),
              size=wx.Size(small_box_dimensions), style=0, value=u'')
        self.input_ftc_window_sec.SetToolTip(
              'Time limit in seconds from the start of the cropped video. '
              'Blank = the whole cropped video.')
        self.text_flies_per_vial = wx.StaticText(id=wx.ID_ANY,
              label=u'Flies / vial:', name='text_flies_per_vial', parent=self.panel1,
              pos=wx.Point(390, row7 + 3), size=wx.Size(65, 22), style=0)
        self.input_flies_per_vial = wx.TextCtrl(id=wx.ID_ANY,
              name=u'input_flies_per_vial', parent=self.panel1, pos=wx.Point(458, row7),
              size=wx.Size(90, 22), style=0, value=u'')
        self.input_flies_per_vial.SetToolTip(
              'Flies loaded per vial: one number for all vials, or a comma-separated '
              'list left to right (e.g. 10,10,9). Lets flies the detector cannot see '
              'still count. Blank = count detected flies only. A vials.txt "n =" line '
              'overrides this per folder.')
        self.text_floor_y = wx.StaticText(id=wx.ID_ANY,
              label=u'Vial floor y (px):', name='text_floor_y', parent=self.panel1,
              pos=wx.Point(560, row7 + 3), size=wx.Size(95, 22), style=0)
        self.input_floor_y = wx.TextCtrl(id=wx.ID_ANY,
              name=u'input_floor_y', parent=self.panel1, pos=wx.Point(658, row7),
              size=wx.Size(medium_box_dimensions), style=0, value=u'')
        self.input_floor_y.SetToolTip(
              'Image y-coordinate of the vial floor, as read off the video axes. '
              'Heights are measured from here. Blank = bottom edge of the ROI.')
        self.text_ftc_hint = wx.StaticText(id=wx.ID_ANY,
              label=u'(Test parameters draws the floor and line on the check frame)',
              name='text_ftc_hint', parent=self.panel1,
              pos=wx.Point(715, row7 + 3), size=wx.Size(230, 30), style=0)

        ## Bottom panels
        self.text_video_path = wx.StaticText(id=wxID_video_path,
              label='Video Path', name='text_video_path', parent=self.panel1,
              pos=wx.Point(10, 240), size=wx.Size(930, 22), style=0)
        self.text_video_path.SetBackgroundColour(wx.Colour(241, 241, 241))

        self.button_test_parameters = wx.Button(id=wxID_test_parameters,
              label=u'Test parameters', name=u'button_test_parameters',
              parent=self.panel1, pos=wx.Point(col1,180), size=wx.Size(140, 22),
              style=wx.ALIGN_CENTER)
        self.button_test_parameters.Bind(wx.EVT_BUTTON, self.OnButton_testParButton,
              id=wxID_test_parameters)

        self.button_reload_video = wx.Button(id=wxID_reload_video,
              label=u'Reload video', name=u'button_reload_video', parent=self.panel1,
              pos=wx.Point(col1 + 160*1, 180), size=wx.Size(140, 22), style=wx.ALIGN_CENTER)
        self.button_reload_video.Bind(wx.EVT_BUTTON, self.OnButton_LoadVideo,
              id=wxID_reload_video)
              
        self.button_store_parameters = wx.Button(id=wxID_store_parameters,
              label=u'Save configuration', name=u'button_store_parameters',
              parent=self.panel1, pos=wx.Point(col1 + 160*2, 180), size=wx.Size(140, 22),
              style=wx.ALIGN_CENTER)
        self.button_store_parameters.Bind(wx.EVT_BUTTON, self.OnButton_strParButton,
              id=wxID_store_parameters)

        ## Text box at the bottom
        self.status_bar = wx.StatusBar(id=wxID_status_bar,
              name='status_bar', parent=self, style=0)
        self.initialize_sizers()
        return   

## Create the GUI
def create(parent, file_path):
    '''Begins running the main GUI object
    ----
    Inputs:
      file_path (str): path to the video file
    ----
    Returns:
      main_gui (object): begins running the main_gui object
    '''
    if args.debug: print('main_gui.create')
    if file_path==None:
        file_path='Select a video to begin'
    return main_gui(parent, file_path)

## Print startup message and setup argument parser
def startup():
    '''Function prints top line and runs argument parsing function.
    ----
    Inputs:
      None
    ----
    Returns:
      args (list): list of arguments passed to program
    '''
    border = '#' * 72
    print('\n')
    print(border)
    print('FreeClimber-FNG v%s' % version)
    print('Built on FreeClimber (Spierer et al., 2020)')
    print('Please cite: https://doi.org/10.5281/zenodo.17647589')
    print('Please cite: https://doi.org/10.1242/jeb.229377')
    print(border)

    args = define_argument_parser()
    return args

def define_argument_parser():
    '''Parses out arguments (which are admittedly few for the GUI. Use the help flag '-h' for more information when running the program.
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

    ## Specifying file (video file)
    parser.add_argument('--video_file', 
                        type=str, 
                        required=False, 
                        help="Path to file")

    ## Debug printing
    parser.add_argument('--debug', 
                        required=False, 
                        default=False, 
                        action='store_true',
                        help="Prints for every function and output being run. Not very pretty")

    ## Final parser specifications
    parser.version = version ## Program version
    parser.add_argument('-v','--version', 
                        action='version')
    args = parser.parse_args()
    return args

## Check video argument arguments
def check_args(args):
    '''Confirms video path is valid, and if not hecks arguments passed to parser and kills program if invalid path.
    ----
    Inputs:
      args (list): Arguments passed to parsed
    ----
    Returns:
      args.video_file (str): Video file path if valid, or None if not.
    '''
    if args.debug: print('main_gui.check_args')
    try:
        if os.path.isfile(args.video_file):
            print('--> Video file: ',args.video_file)
        else:
            print('--> No video file specified, select from Browser...')
            args.video_file = None
    except:
        args.video_file = None
    return args.video_file


## Initializing wx text boxes and input fields
[wxID_text_step_1a,wxID_text_step_1b,wxID_text_step_2,wxID_text_step_3,wxID_text_step_4,wxID_text_step_5,
wxID_panel_1, wxID_text_title, 
wxID_browse_video, wxID_reload_video, wxID_store_parameters,
wxID_test_parameters, wxID_display_parameters,
wxID_video_path, 
wxID_frame_rate,
wxID_status_bar, 
wxID_text_minmass,wxID_input_minmass, 
wxID_text_threshold, wxID_input_threshold,
wxID_text_maxsize, wxID_input_maxsize,
wxID_text_diameter, wxID_input_diameter,
wxID_text_ecc, wxID_input_ecc_low, wxID_input_ecc_high,
wxID_text_x, wxID_input_x,
wxID_text_y, wxID_input_y,
wxID_text_h, wxID_input_h,
wxID_text_w, wxID_input_w,
wxID_text_check_frames, wxID_input_check_frame, 
wxID_text_background_frames,wxID_input_blank_0, wxID_input_blank_n,
wxID_text_crop_frames,wxID_input_crop_0,wxID_input_crop_n,
wxID_text_vials, wxID_input_vials,
wxID_text_window, wxID_input_window,
wxID_text_pixel_to_cm,wxID_input_pixel_to_cm,
wxID_check_box_outlier,wxID_outlier_TB,wxID_outlier_LR,
wxID_text_naming_convention,wxID_input_naming_convention,
wxID_text_vial_id_vars,wxID_input_vial_id_vars,
wxID_text_path_project,wxID_input_path_project,
wxID_input_convert_to_cm_sec,wxID_check_box_ROI] = [wx.ID_ANY for item in range(60)] 


## Basic GUI sizes and spacers
col1,col2,col3,col4,col5 = 10,180,320,525,750
small_box_dimensions = 35,22
medium_box_dimensions = 50,22
box_dimensions = 100,22
large_box_dimensions = 165,22

## Datetime
now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

if __name__=='__main__':
    args = startup()
    file_name = check_args(args)
    app = App(0)
    app.MainLoop()
