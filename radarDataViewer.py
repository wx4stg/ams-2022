import wx

class radarDataViewerFrame(wx.Frame):

    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        ppiPanel = wx.Panel(self, pos=(0,0), size=(256,512))
        rhiPanel = wx.Panel(self, pos=(256,0), size=(256,512))
        ppiTitle = wx.StaticText(ppiPanel, label="Plan Position Indicator")
        rhiTitle = wx.StaticText(rhiPanel, label="Range-Height Indicator")
        ppiAutoLayout = wx.BoxSizer(wx.HORIZONTAL)
        ppiAutoLayout.Add(ppiTitle, wx.SizerFlags().Align(wx.RIGHT))
        ppiPanel.SetSizer(ppiAutoLayout)
        rhiAutoLayout = wx.BoxSizer(wx.HORIZONTAL)
        rhiAutoLayout.Add(rhiTitle, wx.SizerFlags().Align(wx.LEFT))
        rhiPanel.SetSizer(rhiAutoLayout)



if __name__ == "__main__":
    app = wx.App()
    frame = radarDataViewerFrame(None, title="Radar Data Viewer")
    frame.Show()
    app.MainLoop()