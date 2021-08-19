from matplotlib import figure
import wx
import planPositionIndicator
import rangeHeightIndicator
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

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


class radarDataViewerFrame(wx.Frame):

    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        ppiAutoLayout = wx.BoxSizer(wx.VERTICAL)
        ppiPanel = wx.Panel(self, size=(788,820))
        ppiTitle = wx.StaticText(ppiPanel, label="Plan Position Indicator", style=wx.ALIGN_CENTER_HORIZONTAL)
        ppiAutoLayout.Add(ppiTitle, 0, wx.ALIGN_CENTER, 0)
        ppiPlotPanel = mplPanel(planPositionIndicator.plot_radar("TAMU_20210409_0205", True), ppiPanel)
        ppiPlotPanel.SetBackgroundColour("light gray")
        ppiAutoLayout.Add(ppiPlotPanel, 0, wx.EXPAND | wx.ALL, 20)
        ppiPanel.SetSizer(ppiAutoLayout)

        rhiAutoLayout = wx.BoxSizer(wx.VERTICAL)
        rhiPanel = wx.Panel(self, size=(788, 820))
        rhiTitle = wx.StaticText(rhiPanel, label="Range-Height Indicator", style=wx.ALIGN_CENTER_HORIZONTAL)
        rhiAutoLayout.Add(rhiTitle, 0, wx.ALIGN_CENTER, 0)
        rhiPlotPanel = mplPanel(rangeHeightIndicator.plot_crosssection("TAMU_20210409_0205", 350, True), rhiPanel)
        rhiPlotPanel.SetBackgroundColour("light gray")
        rhiAutoLayout.Add(rhiPlotPanel, 0, wx.EXPAND | wx.ALL, 20)
        rhiPanel.SetSizer(rhiAutoLayout)

        twoColumnSizer = wx.BoxSizer(wx.HORIZONTAL)
        twoColumnSizer.Add(ppiPanel, 1, 0, 0)
        twoColumnSizer.Add(rhiPanel, 1, 0, 0)

        

        self.SetSizer(twoColumnSizer)

        self.Fit()



if __name__ == "__main__":
    app = wx.App()
    frame = radarDataViewerFrame(None, title="Radar Data Viewer")
    frame.Show()
    app.MainLoop()