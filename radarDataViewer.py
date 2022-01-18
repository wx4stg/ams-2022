#!/usr/bin/env python3
# wxPython-based GUI for researching plan view vs cross sectional radar data

from os import listdir, getcwd, path
import wx
from wx.core import EVT_BUTTON, EVT_CHECKBOX, EVT_CHOICE, EVT_SLIDER
import planPositionIndicator
import rangeHeightIndicator
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

class mplPanel(wx.Panel):
    def __init__(self, figureToPlot=None, *args, **kw):
        super().__init__(*args, **kw)
        self.SetSize(256, 256)
        if figureToPlot == None:
            errTextAutoLayout = wx.BoxSizer()
            errText = wx.StaticText(self, label="No plot to display.", style=wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL)
            errTextAutoLayout.Add(errText, 1, wx.ALIGN_CENTER, 0)
            self.SetSizer(errTextAutoLayout)
        else:
            self.figure = figureToPlot
            self.canvas = FigureCanvas(self, -1, self.figure)
            self.Fit()
    def updatePlot(self, figureToPlot=None):
        if figureToPlot == None:
            errTextAutoLayout = wx.BoxSizer()
            errText = wx.StaticText(self, label="No plot to display.", style=wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL)
            errTextAutoLayout.Add(errText, 1, wx.ALIGN_CENTER, 0)
            self.SetSizer(errTextAutoLayout)
        else:
            self.figure = figureToPlot
            self.canvas = FigureCanvas(self, -1, self.figure)
            self.Fit()


class radarDataViewerFrame(wx.Frame):

    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self.requestedAz = None
        ppiAutoLayout = wx.BoxSizer(wx.VERTICAL)
        ppiPanel = wx.Panel(self, size=(788,820))
        ppiTitle = wx.StaticText(ppiPanel, label="Plan Position Indicator", style=wx.ALIGN_CENTER_HORIZONTAL)
        ppiAutoLayout.Add(ppiTitle, 0, wx.ALIGN_CENTER, 0)
        ppiPlotPanel = mplPanel(None, ppiPanel, size=(768,768))
        ppiPlotPanel.SetBackgroundColour("light gray")
        self.ppiPlotPanel = ppiPlotPanel
        ppiAutoLayout.Add(ppiPlotPanel, 0, 0, 20)
        ppiPanel.SetSizer(ppiAutoLayout)

        rhiAutoLayout = wx.BoxSizer(wx.VERTICAL)
        rhiPanel = wx.Panel(self, size=(788, 820))
        rhiTitle = wx.StaticText(rhiPanel, label="Range-Height Indicator", style=wx.ALIGN_CENTER_HORIZONTAL)
        rhiAutoLayout.Add(rhiTitle, 0, wx.ALIGN_CENTER, 0)
        rhiPlotPanel = mplPanel(None, rhiPanel, size=(768,768))
        rhiPlotPanel.SetBackgroundColour("light gray")
        self.rhiPlotPanel = rhiPlotPanel
        rhiAutoLayout.Add(rhiPlotPanel, 0, 0, 20)
        rhiPanel.SetSizer(rhiAutoLayout)

        twoColumnSizer = wx.BoxSizer(wx.HORIZONTAL)
        twoColumnSizer.Add(ppiPanel, 1, 0, 0)
        twoColumnSizer.Add(rhiPanel, 1, 0, 0)

        controlsPanelSizer = wx.BoxSizer(wx.HORIZONTAL)
        controlsPanel = wx.Panel(self, size=(788, 100))
        radarDataDir = path.join(getcwd(), "radarData") 
        availTimesDropDown = wx.Choice(controlsPanel, choices=sorted(listdir(radarDataDir)), style=wx.ALIGN_CENTER_VERTICAL)
        availTimesDropDown.Bind(EVT_CHOICE, self.onTimeSel)
        controlsPanelSizer.Add(availTimesDropDown, 1, wx.EXPAND | wx.ALL, 20)
        shouldShowLtgBox = wx.CheckBox(controlsPanel, label="Show lightning?", style=wx.ALIGN_CENTER_VERTICAL)
        shouldShowLtgBox.Bind(EVT_CHECKBOX, self.onSetShowLtg)
        self.shouldPlotLightning = False
        controlsPanelSizer.Add(shouldShowLtgBox, 1, wx.EXPAND | wx.ALL, 20)
        shouldShowRptsBox = wx.CheckBox(controlsPanel, label="Show Hail Reports?", style=wx.ALIGN_CENTER_VERTICAL)
        shouldShowRptsBox.Bind(EVT_CHECKBOX, self.onSetShowReports)
        self.shouldPlotReports = False
        controlsPanelSizer.Add(shouldShowRptsBox, 1, wx.EXPAND | wx.ALL, 20)
        radiusText = wx.StaticText(controlsPanel, label="PPI map radius (km):", style=wx.ALIGN_CENTER_HORIZONTAL)
        controlsPanelSizer.Add(radiusText, 1, wx.EXPAND | wx.ALL, 20)
        radSlider = wx.Slider(controlsPanel, value=30, minValue=10, maxValue=300, style= wx.SL_HORIZONTAL | wx.SL_LABELS, size=(100, 100))
        self.requestedRad = 30
        radSlider.Bind(EVT_SLIDER, self.onRadSliderScroll)
        controlsPanelSizer.Add(radSlider, 10, wx.EXPAND | wx.ALL, 20)
        rrText = wx.StaticText(controlsPanel, label="Range ring step (km):", style=wx.ALIGN_CENTER_HORIZONTAL)
        controlsPanelSizer.Add(rrText, 1, wx.EXPAND | wx.ALL)
        rrSlider = wx.Slider(controlsPanel, value=10, minValue=0, maxValue=300, style=wx.SL_HORIZONTAL | wx.SL_LABELS, size=(100, 100))
        self.requestedRRStep = 10
        rrSlider.Bind(EVT_SLIDER, self.onrrSliderScroll)
        controlsPanelSizer.Add(rrSlider, 10, wx.EXPAND | wx.ALL, 20)
        azText = wx.StaticText(controlsPanel, label="Lat for cross-section:", style=wx.ALIGN_CENTER_HORIZONTAL)
        controlsPanelSizer.Add(azText, 1, wx.EXPAND | wx.ALL, 20)
        azSlider = wx.Slider(controlsPanel, value=3025, minValue=3025, maxValue=3225, style= wx.SL_HORIZONTAL | wx.SL_LABELS, size=(100, 100))
        azSlider.Bind(EVT_SLIDER, self.onAzSliderScroll)
        controlsPanelSizer.Add(azSlider, 10, wx.EXPAND | wx.ALL, 20)
        redrawBtn = wx.Button(controlsPanel, label="Redraw", size=(100, 25))
        redrawBtn.Bind(EVT_BUTTON, self.onRedrawBtnPress)
        controlsPanelSizer.Add(redrawBtn, 1, wx.EXPAND | wx.ALL, 20)
        exportBtn = wx.Button(controlsPanel, label="Export", size=(100, 25))
        exportBtn.Bind(EVT_BUTTON, self.onExportBtnPress)
        controlsPanelSizer.Add(exportBtn, 1, wx.EXPAND | wx.ALL, 20)
        controlsPanel.SetSizer(controlsPanelSizer)
        
        plotsAndControlsSizer = wx.BoxSizer(wx.VERTICAL)
        plotsAndControlsSizer.Add(twoColumnSizer, 0, wx.ALIGN_CENTER, 10)
        plotsAndControlsSizer.Add(controlsPanel, 0, wx.EXPAND | wx.BOTTOM, 1)
        self.SetSizer(plotsAndControlsSizer)
        self.Fit()

    def onTimeSel(self, event):
        obj = event.GetEventObject()
        self.requestedFile = obj.GetString(obj.GetSelection())
        self.drawPPI()
    def onSetShowLtg(self, event):
        obj = event.GetEventObject()
        self.shouldPlotLightning = obj.GetValue()
    def onSetShowReports(self, event):
        obj = event.GetEventObject()
        self.shouldPlotReports = obj.GetValue()
    def onAzSliderScroll(self, event):
        self.requestedAz = event.GetEventObject().GetValue()
    def onRadSliderScroll(self, event):
        self.requestedRad = event.GetEventObject().GetValue()
    def onrrSliderScroll(self, event):
        self.requestedRRStep = event.GetEventObject().GetValue()
    def onRedrawBtnPress(self, event):
        self.drawPPI()
        self.drawRHI()
    def onExportBtnPress(self, event):
        planPositionIndicator.plot_radar(self.requestedFile, False, 30, None, str("exports/ppi_"+self.requestedFile+".png"))
        rangeHeightIndicator.plot_crosssection(self.requestedFile, str("exports/rhi_"+self.requestedFile+".png"), self.requestedAz/100, False, 200, 10, self.shouldPlotLightning)
    def drawPPI(self):
        if self.requestedFile is not None:
            self.ppiPlotPanel.updatePlot(planPositionIndicator.plot_radar(self.requestedFile, None, True, self.requestedRad, self.requestedRRStep, self.requestedAz, self.shouldPlotReports, self.shouldPlotLightning))
    def drawRHI(self):
        if self.requestedFile is not None:
            if self.requestedAz is not None:
                rangeHeightIndicator.plot_crosssection(self.requestedFile, str("exports/rhi_"+self.requestedFile+".png"), self.requestedAz/100, True, 200, 10, self.shouldPlotLightning)

if __name__ == "__main__":
    app = wx.App()
    frame = radarDataViewerFrame(None, title="Radar Data Viewer")
    frame.Show()
    app.MainLoop()
