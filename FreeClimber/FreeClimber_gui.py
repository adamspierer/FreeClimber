#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import wx
import sys
import matplotlib
matplotlib.use('WXAgg')
import DectFrame
import detector


class App(wx.App):
    def OnInit(self):
        '''Create the main window and insert the custom frame'''
        frame = DectFrame.create(None, file_name)
        self.SetTopWindow(frame)
        frame.Show(True)
        return True

if __name__=='__main__':
    if len(sys.argv)>1:
    	file_name = sys.argv[1]
    else:
    	file_name = None
    app = App(0)
    app.MainLoop()
