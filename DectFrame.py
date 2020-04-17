#!/usr/bin/env python
# -*- coding: utf-8 -*- 
#Boa:Frame:DectFrame

import wx
import os
import sys
import time
import numpy as np
import pandas as pd

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from detector import detector
import vial_detector_main

print('')

class DectFrame(wx.Frame):    
    def __init__(self, parent, file_path):
        print('-------->  __init__')
        self._init_controls(parent)
        self.box_sizer.Add(self.panel1, 0, wx.EXPAND)

        self.pressed=False
        ## Initialize GUI plot
        self.figure = Figure()
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.box_sizer.Add(self.canvas, 1, wx.EXPAND)

        ## Initialize bottom text bar
        self.box_sizer.Add(self.status_bar, 0, border=0, flag=0)
        self.status_bar.SetStatusText('Ready', 0)
        rect = self.status_bar.GetFieldRect(0)
    

        ## Set up names
        self.video_file = file_path
        
        ## Check input file is valid
        if not os.path.isfile(file_path):
            openFileDialog = wx.FileDialog(self, "Open Video file", "", "",
                                       "Video files (*.*)|*.*", 
                                       wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

            if openFileDialog.ShowModal() == wx.ID_CANCEL:
                print('\n\nExiting Program. Must select a video from box or enter path into command line\n\n')
                raise SystemExit

            self.video_file = openFileDialog.GetPath()        
            file_path = self.video_file
        ## Passing individual inputs from GUI to a list
        self.update_names()
        self.input_names = ['x','y','w','h','frame_0','frame_n','blank_0','blank_n',
        'threshold','diameter','minmass','maxsize','vials','window','pixel_to_cm','frame_rate',
        'vvars','naming_convention','path_project','file_suffix','convert_to_cm']
        self.parameter_names = ['self.input_'+item for item in self.input_names]
        self.input_values = [self.input_x,self.input_y,self.input_w,self.input_h,
                             self.input_frame_0,self.input_frame_n,self.input_blank_0,
                             self.input_blank_n,self.input_threshold,self.input_diameter,
                             self.input_minmass,self.input_maxsize,self.input_vials,
                             self.input_window,self.input_pixel_to_cm,self.input_frame_rate,
                             self.input_vvars,self.input_naming_convention,self.input_path_project,
                             self.input_file_suffix,self.input_convert_to_cm]

        ## Enable all buttons
        button_list = ['browse_video','reload_video','test_parameters','store_parameters','run_analysis']
        for button in button_list:
            exec('self.button_' + button + '.Enable(True)')

        self.load_video()
        print('End of init')
        return

    def update_variables(self):
        print('---------> update_variables')
        variables = []
        
#         self.input_convert_to_cm = self.input_convert_to_cm.GetValue()
        for item,jtem in zip(self.input_names[:17],self.input_values[:17]):
            phrase = str(item + '='+ jtem.GetValue())
            print('    '+phrase)
            variables.append(phrase)
        for item,jtem in zip(self.input_names[17:19],self.input_values[17:19]):
            phrase = str(item + '="'+ jtem.GetValue()+'"')
            print('    '+phrase)
            variables.append(phrase)
        for item,jtem in zip(self.input_names[19:20],self.input_values[19:20]):
            phrase = str(item + '="'+ str(jtem)+'"')
            print(phrase)
            variables.append(phrase)
        for item,jtem in zip(self.input_names[20:],self.input_values[20:]):
            phrase = str(item + '=%s' % str(jtem.GetValue()))
            print('    '+phrase)
            variables.append(phrase)

        return variables

    def load_video(self):
        print('--------> load_video')
        #Load the video, set Cursor to be busy
        self.figure.clear()
        self.axes   = [self.figure.add_subplot(111), ]
        
        ## Confirm file is a path, or folder has specified suffix
        status = self.check_specified_video()
        if status:
            print("Loading:", self.video_file)
#             self.input_frame_rate.SetValue('25')
            self.update_names()
            for item in self.input_values[:4]:
                item.SetEditable(True)
                item.Enable(True)
            self.checkBox_fixed_ROI.Enable(True)
            self.input_convert_to_cm.Enable(True)
            wx.BeginBusyCursor()
            try:
                variables = self.update_variables()
                self.detector = detector(self.video_file,
                                        gui=True,
                                        variables = variables)
#                 self.input_frame_rate.SetValue(str(self.detector.frame_rate))
            
                self.axes[0].imshow(self.detector.image_stack[0])
                self.figure.canvas.draw()
            finally:
                wx.EndBusyCursor()

            self.canvas.Bind(wx.EVT_ENTER_WINDOW, self.ChangeCursor)
            self.canvas.mpl_connect('button_press_event', self.draw_rectangle)
            self.canvas.mpl_connect('button_release_event', self.on_release)
            self.canvas.mpl_connect('motion_notify_event', self.on_motion)
            self.rect = Rectangle((0,0), 1, 1, fill=False, ec='r')
            self.axes[0].add_patch(self.rect)


            ## Auto-set GUI parameters form video
#             self.input_frame_rate.SetValue(str(self.input_frame_rate))
            self.input_blank_0.SetValue('0')
            self.input_blank_n.SetValue(str(self.detector.n_frames))
            self.input_frame_0.SetValue('0')
            self.input_path_project.SetValue(self.folder)
            self.input_naming_convention.SetValue(self.name)
            self.input_vvars.SetValue(str(len(self.input_naming_convention.GetValue().split('_'))))
            ## Display the 0th and frame corresponding with (most likely) t = 2 seconds
            try:
                self.input_frame_rate = int(self.input_frame_rate.GetValue())
            except:
                pass
            if self.detector.n_frames < self.input_frame_rate*2:            
                self.input_frame_n.SetValue(str(self.detector.n_frames))
            else:
                self.input_frame_n.SetValue(str(self.input_frame_rate*2))
            
            ## Try to make the local linear regression window size 2 seconds, but if not then 35% of the frames in the video
            if self.detector.n_frames < self.input_frame_rate*2:
                self.input_window.SetValue(str(int(len(self.detector.image_stack) * .35)))                
            else:
                self.input_window.SetValue(str(int(self.input_frame_rate)*2))

            ## Enable Test parameter button if disabled from prior testing
            self.button_test_parameters.Enable(True)
        
            self.x0 = 0
            self.y0 = 0
            self.x1 = self.detector.width
            self.y1 = self.detector.height
        
            self.update_ROIdisp()
            self.canvas.draw()
        else:
            return

    def update_names(self):
        print('--------> update_names')
        self.text_video_path.SetLabelText(self.video_file)
        self.folder, self.name        = os.path.split(self.video_file)
        print(self.name)
        self.name,   self.input_file_suffix = self.name.split('.')

        self.name_noext  = os.path.join(self.folder,self.name)
        self.path_data   = self.name_noext + '.raw.csv'
        self.path_filter = self.name_noext + '.filter.csv'
        self.path_plot   = self.name_noext + '.diag.png'
        self.path_slope  = self.name_noext + '.slopes.csv'
        
        if self.input_path_project == '':
            self.input_path_project = self.folder


    def check_specified_video(self):
        print('--------> check_specified_video')
        if os.path.isfile(self.video_file):
            self.button_reload_video.Enable(True)
            self.update_names()
            
#             self.button_reload_video.Enable(False)
            self.button_test_parameters.Enable(True)
            
            self.input_file_suffix = '.' + self.video_file.split('/')[-1].split('.')[-1]
            return True

        elif file_path == None:
            self.button_browse_video.Enable(True)
            return False

        else:
            self.video_file = "Invalid file. Please change the file path"
            self.button_browse_video.Enable(True)
            return False


    def ChangeCursor(self, event):
        '''Change cursor into crosshair type when enter the plot area'''
        self.canvas.SetCursor(wx.Cursor(wx.CURSOR_CROSS))


    def draw_rectangle(self, event):
        '''Draw ROI rectangle'''
        self.pressed = True
        if self.checkBox_fixed_ROI.Enabled:
            try:
                self.x0 = int(event.xdata)
                self.y0 = int(event.ydata)
                if self.checkBox_fixed_ROI.GetValue():
                    self.x1 = self.x0 + int(eval(self.input_w.GetValue()))
                    self.y1 = self.y0 + int(eval(self.input_h.GetValue()))
                    self.rect.set_width(self.x1 - self.x0)
                    self.rect.set_height(self.y1 - self.y0)
                    self.rect.set_xy((self.x0, self.y0))
                    self.canvas.draw()
                self.input_x.SetValue(str(self.x0))
                self.input_y.SetValue(str(self.y0))
                self.input_h.SetValue(str(self.rect.get_height()))
                self.input_w.SetValue(str(self.rect.get_width()))
            except:
                pass


    def on_release(self, event):
        '''When mouse is on plot and button is released, redraw ROI rectangle, update ROI values'''
        self.pressed = False
        if self.checkBox_fixed_ROI.Enabled:
            if self.checkBox_fixed_ROI.GetValue():
                pass
            else:
                self.redraw_rect(event)
            self.update_ROIdisp()


    def on_motion(self, event):
        '''If the mouse is on plot and if the mouse button is pressed, redraw ROI rectangle'''
        if self.pressed & self.checkBox_fixed_ROI.Enabled & (not self.checkBox_fixed_ROI.GetValue()):
            # redraw the rect
            self.redraw_rect(event)
            self.update_ROIdisp()


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
        except:
            pass


    def update_ROIdisp(self):
        '''Update ROI display in UI'''
        self.input_x.SetValue(str(self.x0))
        self.input_y.SetValue(str(self.y0))
        self.input_h.SetValue(str(int(self.y1) - int(self.y0)))
        self.input_w.SetValue(str(int(self.x1) - int(self.x0)))
       

    def OnButton_runAnaButton(self, event):
        '''Analysis the whole video, output the result pandas.dataframe and analysis config file
        in the same folder (as the video file).
        Plot diagnostic plots.
        '''
        print('-------->  OnButton_runAnaButton')
        variables = self.update_variables()
        time_begin = time.time()
        self.status_bar.SetStatusText("Running analysis...",0)
        self.button_run_analysis.Enable(False)
        self.save_parameter()

        self.figure.clear()
        self.axes = self.figure.add_subplot(111)
        wx.BeginBusyCursor()

        ## Running the main program
        try:  
#             free_climber_main.free_climber.process(File = self.video_file,variables=variables)
#             raise SystemExit
            os.system('python ./vial_detector_main.py %s False' % (self.path_parameters)) ## HELP
        finally:
            self.status_bar.SetStatusText("Analysis complete. Reload video to test new parameters, try a different video, or batch process all videos using command line",0)
            wx.EndBusyCursor()

        ## Plotting the resulting diagnostic plot
        ##   --> Would rather this be plot directly into the GUI, rather than loading
        ##       a picture to the GUI that has poor resolution
        image = plt.imread(self.path_plot)
        self.axes.imshow(image)
        self.axes.spines['top'].set_visible(False)
        self.axes.spines['right'].set_visible(False)
        self.axes.xaxis.set_visible(False)
        self.axes.yaxis.set_visible(False)
        self.figure.canvas.draw()       

        self.button_reload_video.Enable(True)
        self.button_store_parameters.Enable(True)
        self.button_run_analysis.Enable(True)
        self.button_test_parameters.Enable(True)
        return

    def OnButton_testParButton(self, event):
        '''plot the decection result and compare two frames, an early and a late Frame
        So the parameter can be fine tuned.'''
        print('-------->  OnButton_testParButton')

        #Prep the parameters
        variables = self.update_variables()
        self.checkBox_fixed_ROI.Enable(False)

        ## Set up figure for plots
        self.figure.clear()
        self.axes = [self.figure.add_subplot(231),
                     self.figure.add_subplot(232),
                     self.figure.add_subplot(233),
                     self.figure.add_subplot(234),
                     self.figure.add_subplot(235),
                     self.figure.add_subplot(236)]

        wx.BeginBusyCursor()
        try:
            self.detector.parameter_testing(variables, self.axes)
        finally:
            wx.EndBusyCursor()
        
        self.figure.tight_layout()
        self.figure.canvas.draw()

        #Once the ROI is defined, disable further changes
        self.button_reload_video.Enable(True)
        self.button_store_parameters.Enable(True)
        self.button_run_analysis.Enable(True)
        self.button_test_parameters.Enable(False)


    def OnButton_strParButton(self, event):
        print('-------->  OnButton_strParButton')
        self.save_parameter()

    def set_config_file(self):
        '''DocString'''
        ## Figure out where to save configuration file
        if os.path.isdir(self.input_path_project):
            if not self.input_path_project.endswith('/'):
                self.input_path_project = self.input_path_project + '/'
            self.path_parameters = self.input_path_project + self.name + '.cfg'
        else:
            self.path_parameters = self.path_noext+'.cfg'
        return self.path_parameters

    def save_parameter(self):
        print('-------->  save_parameter')
        '''
        Save parameters as python list.
        New parameter sets appended to the config file.
        Each parameter sets come with a comment line, contain the datetime of analysis
        '''

        variables = self.update_variables()
        try: self.input_path_project = self.input_path_project.GetValue()
        except: pass

        self.path_parameters = self.set_config_file()

        ## Printing output to configuration file
        print('Saving parameters to:', self.path_parameters)
        with open(self.path_parameters, 'w') as f:
            print('## Vial Detector ##', file=f)
        f.close()
        with open(self.path_parameters,'a') as f:
            print('## Generated from file:'+self.video_file,file = f)
            print('##     @ ' + str(time.ctime()), file = f)
            print('##',file = f)
            print('## Analysis parameters:',file = f)
            for item in variables:
                print(item,file = f)
        f.close()
        print("Configuration settings saved")
        return

    def OnButton_Browse(self, event):
        print('-------->  OnButton_Browse')
#         self.input_frame_rate.SetValue('1')
        openFileDialog = wx.FileDialog(self, "Open Video file", "", "",
                                       "Video files (*.*)|*.*", 
                                       wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        if openFileDialog.ShowModal() == wx.ID_CANCEL:
            pass
        else:
            self.video_file = openFileDialog.GetPath()
            self.update_names()
            self.load_video()
        self.figure.clear()
        return


    def OnButton_LoadVideo(self,event):
        print('-------->  OnButton_LoadVideo')
        self.load_video()
        return

        
    def _init_sizers(self):
        print('-------->  _init_sizers')
        # Generated method, do not edit
        self.box_sizer = wx.BoxSizer(orient=wx.VERTICAL)
        self.SetSizer(self.box_sizer)

    def _init_controls(self, prnt):
        print('-------->  _init_controls')
        # Generated method, do not edit
        wx.Frame.__init__(self, id=wxID_text_title, name='', parent=prnt,
              pos=wx.Point(100, 30), size=wx.Size(950, 759),
              style=wx.DEFAULT_FRAME_STYLE,
              title='Background adjusted particle detection')
        self.SetClientSize(wx.Size(950, 737))

        ######
        ## Inputs for ROI Rectangle
        self.panel1 = wx.Panel(id=wxID_panel_1, name='panel1',
              parent=self, pos=wx.Point(0, 0), size=wx.Size(950, 231),
              style=wx.TAB_TRAVERSAL)
       
        ## Step 1 boxes
        self.text_step_1a = wx.StaticText(id=wxID_text_step_1a,
              label=u'Step 1a: Find video', name='text_step_1a', parent=self.panel1, 
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
              parent=self.panel1, pos=wx.Point(col1, 115), size=wx.Size(medium_box_dimensions),
              style=0)
        self.input_pixel_to_cm = wx.TextCtrl(id=wxID_input_pixel_to_cm,
              name=u'input_pixel_to_cm', parent=self.panel1, pos=wx.Point(col1 + 95, 115), 
              size=wx.Size(medium_box_dimensions), style=0, value=u"1")
        
        ## Frame Rate    
        self.text_frame_rate = wx.StaticText(id=wxID_frame_rate,
              label=u'Frames / sec:', name='text_frame_rate',
              parent=self.panel1, pos=wx.Point(col1,85), size=wx.Size(box_dimensions),
              style=0)
        self.input_frame_rate = wx.TextCtrl(id=wxID_frame_rate,
              name=u'input_frame_rate', parent=self.panel1, pos=wx.Point(col1 + 95,85),
               size=wx.Size(medium_box_dimensions), style=0, value='25')
        
        ## Check box to convert final slope to cm
        self.input_convert_to_cm = wx.CheckBox(id=wxID_input_convert_to_cm,
              label=u'Convert to cm', name=u'input_convert_to_cm',
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
              parent=self.panel1, pos=wx.Point(col2, 130), size=wx.Size(250,22),
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
              name=u'input_diameter', parent=self.panel1, pos=wx.Point(col3 + 90,30), 
              size=wx.Size(medium_box_dimensions), style=0, value=u'7')
        
        ## Maximum spot diameter
        self.text_maxsize = wx.StaticText(id=wxID_text_maxsize,
              label=u'MaxDiameter:', name='text_maxsize', parent=self.panel1,
              pos=wx.Point(col3,55), size=wx.Size(medium_box_dimensions), style=0)
        self.input_maxsize = wx.TextCtrl(id=wxID_input_maxsize,
              name=u'input_maxsize', parent=self.panel1, pos=wx.Point(col3 + 90,55), 
              size=wx.Size(medium_box_dimensions), style=0, value=u'11')

        ## Minimum spot 'mass'
        self.text_threshold = wx.StaticText(id=wxID_text_minmass,
              label=u'MinMass:', name='text_threshold', parent=self.panel1,
              pos=wx.Point(col3,80), size=wx.Size(medium_box_dimensions), style=0)
        self.input_minmass = wx.TextCtrl(id=wxID_input_minmass,
              name=u'input_minmass', parent=self.panel1, pos=wx.Point(col3 + 90,80), 
              size=wx.Size(medium_box_dimensions), style=0, value=u'100')

        ## Spot threshold
        self.text_threshold = wx.StaticText(id=wxID_text_threshold,
              label=u'Threshold:', name='text_threshold', parent=self.panel1,
              pos=wx.Point(col3, 105), size=wx.Size(medium_box_dimensions), style=0)
        self.input_threshold = wx.TextCtrl(id=wxID_input_threshold,
              name=u'input_threshold', parent=self.panel1, pos=wx.Point(col3 + 90,105),
              size=wx.Size(medium_box_dimensions), style=0, value=u'"auto"')
        
        ## Check frames
        self.text_step_4 = wx.StaticText(id=wxID_text_step_4,
              label=u'Step 4: Additional parameters', name='text_step_4', parent=self.panel1, 
              pos=wx.Point(col4,10),
              size=wx.Size(100,22),style = wx.ALIGN_LEFT)

        self.text_check_frames = wx.StaticText(id=wxID_text_check_frames,
              label=u'Check Frames:', name='text_check_frames', parent=self.panel1,
              pos=wx.Point(col4, 55), size=wx.Size(115, 17), style=0)

        self.input_frame_0 = wx.TextCtrl(id=input_frame_0,
              name=u'input_frame_0', parent=self.panel1, pos=wx.Point(col4 + 130, 55),
               size=wx.Size(small_box_dimensions), style=0, value=u'')

        self.input_frame_n = wx.TextCtrl(id=wxID_input_frame_n,
              name=u'input_frame_n', parent=self.panel1, pos=wx.Point(col4 + 170, 55), 
              size=wx.Size(small_box_dimensions), style=0, value=u'')
        
        ## Background frames
        self.text_background_frames = wx.StaticText(id=wxID_text_background_frames,
              label=u'Background Frames:', name='text_background_frames',
              parent=self.panel1, pos=wx.Point(col4, 30), size=wx.Size(medium_box_dimensions),
              style=0)
        self.input_blank_0 = wx.TextCtrl(id=wxID_input_blank_0,
              name=u'input_blank_0', parent=self.panel1, pos=wx.Point(col4+ 130, 30),
            size=wx.Size(small_box_dimensions), style=0, value=u'')

        self.input_blank_n = wx.TextCtrl(id=wxID_input_blank_n,
              name=u'input_blank_n', parent=self.panel1, pos=wx.Point(col4 + 170,30),
              size=wx.Size(small_box_dimensions), style=0, value=u'')
        
        ## Vials
        self.text_vials = wx.StaticText(id=wxID_text_vials,
              label=u'Number of vials:', name='text_vials',
              parent=self.panel1, pos=wx.Point(col4, 80), size=wx.Size(133, 22),
              style=0)
        self.input_vials = wx.TextCtrl(id=wxID_input_vials,
              name=u'input_vials', parent=self.panel1, pos=wx.Point(col4 + 130, 80), 
              size=wx.Size(small_box_dimensions), style=0, value=u'1')
        
        ## Window size
        self.text_window = wx.StaticText(id=wxID_text_window,
              label=u'Window size:', name='text_window',
              parent=self.panel1, pos=wx.Point(col4, 105), size=wx.Size(133, 22),
              style=0)
        self.input_window = wx.TextCtrl(id=wxID_input_window,
              name=u'input_window', parent=self.panel1, pos=wx.Point(col4 + 130, 105), 
              size=wx.Size(small_box_dimensions), style=0, value='1')    
          

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
              name=u'input_naming_convention', parent=self.panel1, pos=wx.Point(col5 + 20, 50), 
              size=wx.Size(large_box_dimensions), style=0, value='') 
        
        ## Variables
        self.text_vvars = wx.StaticText(id=wxID_text_vvars,
              label=u'Video vars:', name='text_vvars',
              parent=self.panel1, pos=wx.Point(col5, 80), size=wx.Size(133, 22),
              style=0)
        self.input_vvars = wx.TextCtrl(id=wxID_input_vvars,
              name=u'input_vvars', parent=self.panel1, pos=wx.Point(col5 + 130, 80), 
              size=wx.Size(small_box_dimensions), style=0, value=u'2')  
        
        self.text_path_project = wx.StaticText(id=wxID_text_path_project,
              label=u"Project path:", name='text_path_project',
              parent=self.panel1, pos=wx.Point(col5, 105), size=wx.Size(medium_box_dimensions),
              style=0)
        self.input_path_project = wx.TextCtrl(id=wxID_input_path_project,
              name=u'input_path_project', parent=self.panel1, pos=wx.Point(col5, 125), 
              size=wx.Size(large_box_dimensions), style=0, value='') 


        ## Bottom panels
        self.text_video_path = wx.StaticText(id=wxID_video_path,
              label='Video Path', name='text_video_path', parent=self.panel1,
              pos=wx.Point(10, 205), size=wx.Size(930, 22), style=0)
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

        self.button_run_analysis = wx.Button(id=wxID_run_analysis,
              label=u'Run analysis', name=u'button_run_analysis', parent=self.panel1,
              pos=wx.Point(col1 + 160*3, 180), size=wx.Size(140, 22), style=wx.ALIGN_CENTER)
        self.button_run_analysis.Bind(wx.EVT_BUTTON, self.OnButton_runAnaButton,
              id=wxID_run_analysis)

#         self.button_display_parameters = wx.Button(id=wxID_display_parameters,
#               label=u'Display Parameters', name=u'button_display_parameters', parent=self.panel1,
#               pos=wx.Point(col1 + 160*4, 180), size=wx.Size(140, 22), style=wx.ALIGN_CENTER)
#         self.button_display_parameters.Bind(wx.EVT_BUTTON, self.OnButton_display_parameters,
#               id=wxID_display_parameters)
              
        self.status_bar = wx.StatusBar(id=wxID_status_bar,
              name='status_bar', parent=self, style=0)
        self._init_sizers()


[wxID_text_step_1a,wxID_text_step_1b,wxID_text_step_2,wxID_text_step_3,wxID_text_step_4,wxID_text_step_5,
wxID_panel_1, wxID_text_title, 
wxID_browse_video, wxID_reload_video,
wxID_run_analysis, wxID_store_parameters,
wxID_test_parameters, wxID_display_parameters,
wxID_video_path, 
wxID_frame_rate,
wxID_status_bar, 
wxID_text_minmass,wxID_input_minmass, 
wxID_text_threshold, wxID_input_threshold,
wxID_text_maxsize, wxID_input_maxsize,
wxID_text_diameter, wxID_input_diameter,
wxID_text_x, wxID_input_x,
wxID_text_y, wxID_input_y,
wxID_text_h, wxID_input_h,
wxID_text_w, wxID_input_w,
wxID_text_check_frames, input_frame_0, wxID_input_frame_n,
wxID_text_background_frames,wxID_input_blank_0, wxID_input_blank_n,
wxID_text_vials, wxID_input_vials,
wxID_text_window, wxID_input_window,
wxID_text_pixel_to_cm,wxID_input_pixel_to_cm,
wxID_text_naming_convention,wxID_input_naming_convention,
wxID_text_vvars,wxID_input_vvars,
wxID_text_path_project,wxID_input_path_project,
wxID_input_convert_to_cm,wxID_check_box_ROI,
wxID_vert_line_1] = [wx.NewId() for _init_controls in range(54)]


## Basic GUI sizes and spacers
col1,col2,col3,col4,col5 = 10,180,320,525,750

small_box_dimensions = 35,22
medium_box_dimensions = 50,22
box_dimensions = 100,22
large_box_dimensions = 145,22

def create(parent, file_path):
    print('-------->  create')
    if file_path==None:
        file_path='Select a video to begin'
    return DectFrame(parent, file_path)
