import wx

class mplPanel(wx.Panel):
    def __init__(self, figureToPlot=None, *args, **kw):
        super().__init__(*args, **kw)
        if figureToPlot == None:
            errTextAutoLayout = wx.BoxSizer()
            errText = wx.StaticText(self, label="No plot to display.", style=wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL)
            errTextAutoLayout.Add(errText, wx.SizerFlags().Align(wx.CENTER))
            self.SetSizer(errTextAutoLayout)


class radarDataViewerFrame(wx.Frame):

    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self.SetSize(512,512)
        ppiPanel = wx.Panel(self, pos=(0,0), size=(256,512))
        ppiTitle = wx.StaticText(ppiPanel, label="Plan Position Indicator", style=wx.ALIGN_CENTER_HORIZONTAL)
        ppiTitle.SetSize(256, ppiTitle.GetClientSize().GetHeight())
        ppiAutoLayout = wx.BoxSizer(wx.VERTICAL)
        ppiAutoLayout.Add(ppiTitle, wx.SizerFlags().Align(wx.CENTER))
        ppiPlot = mplPanel(None, self)
        ppiPlot.SetBackgroundColour("light gray")
        ppiAutoLayout.Add(ppiPlot)
        ppiPanel.SetSizer(ppiAutoLayout)

        rhiPanel = wx.Panel(self, pos=(256,0), size=(256,512))
        rhiTitle = wx.StaticText(rhiPanel, label="Range-Height Indicator", style=wx.ALIGN_CENTER_HORIZONTAL)
        rhiTitle.SetSize(256, rhiTitle.GetClientSize().GetHeight())
        rhiAutoLayout = wx.BoxSizer(wx.VERTICAL)
        rhiAutoLayout.Add(rhiTitle, wx.SizerFlags().Align(wx.CENTER))
        rhiPanel.SetSizer(rhiAutoLayout)



if __name__ == "__main__":
    app = wx.App()
    frame = radarDataViewerFrame(None, title="Radar Data Viewer")
    frame.Show()
    app.MainLoop()