#!/usr/bin/env python3
# Next-gen HDWX radar plotting script
# Created 7 July 2021 by Sam Gardner <stgardner4@tamu.edu>

from datetime import datetime as dt
from datetime import timedelta, timezone
import xarray as xr
import pyart
from matplotlib import pyplot as plt
from os import path, getcwd, listdir
from cartopy import crs as ccrs
from metpy.plots import ctables
from metpy.plots import USCOUNTIES
import numpy as np
import warnings
import multiprocessing as mp
from matplotlib import image as mpimage
import pandas as pd
from pyxlma.lmalib.io import read as readlma


def plot_ppi_map_modified(
            rmd, field, sweep=0, mask_tuple=None,
            vmin=None, vmax=None, cmap=None, norm=None, mask_outside=False,
            title=None, title_flag=True,
            colorbar_flag=True, colorbar_label=None, ax=None, fig=None,
            lat_lines=None, lon_lines=None, projection=None,
            min_lon=None, max_lon=None, min_lat=None, max_lat=None,
            width=None, height=None, lon_0=None, lat_0=None,
            resolution='110m', shapefile=None, shapefile_kwargs=None,
            edges=True, gatefilter=None,
            filter_transitions=True, embelish=True, raster=False,
            ticks=None, ticklabs=None, alpha=None):
        # parse parameters
        ax, fig = pyart.graph.common.parse_ax_fig(ax, fig)
        vmin, vmax = pyart.graph.common.parse_vmin_vmax(rmd._radar, field, vmin, vmax)
        cmap = pyart.graph.common.parse_cmap(cmap, field)
        if lat_lines is None:
            lat_lines = np.arange(30, 46, 1)
        if lon_lines is None:
            lon_lines = np.arange(-110, -75, 1)
        lat_0 = rmd.loc[0]
        lon_0 = rmd.loc[1]

        # get data for the plot
        data = rmd._get_data(
            field, sweep, mask_tuple, filter_transitions, gatefilter)
        x, y = rmd._get_x_y(sweep, edges, filter_transitions)

        # mask the data where outside the limits
        if mask_outside:
            data = np.ma.masked_outside(data, vmin, vmax)

        # initialize instance of GeoAxes if not provided
        if hasattr(ax, 'projection'):
            projection = ax.projection
        else:
            if projection is None:
                # set map projection to LambertConformal if none is specified
                projection = ccrs.LambertConformal(
                    central_longitude=lon_0, central_latitude=lat_0)
                warnings.warn("No projection was defined for the axes."
                              + " Overridding defined axes and using default "
                              + "axes.", UserWarning)
            ax = plt.axes(projection=projection)

        if min_lon:
            ax.set_extent([min_lon, max_lon, min_lat, max_lat],
                          crs=ccrs.PlateCarree())
        elif width:
            ax.set_extent([-width/2., width/2., -height/2., height/2.],
                          crs=rmd.grid_projection)

        # plot the data
        if norm is not None:  # if norm is set do not override with vmin/vmax
            vmin = vmax = None
        pm = ax.pcolormesh(x * 1000., y * 1000., data, alpha=alpha,
                           vmin=vmin, vmax=vmax, cmap=cmap,
                           norm=norm, transform=rmd.grid_projection)

        # plot as raster in vector graphics files
        if raster:
            pm.set_rasterized(True)

        if title_flag:
            rmd._set_title(field, sweep, title, ax)

        # add plot and field to lists
        rmd.plots.append(pm)
        rmd.plot_vars.append(field)

        if colorbar_flag:
            rmd.plot_colorbar(
                mappable=pm, label=colorbar_label, field=field, fig=fig,
                ax=ax, ticks=ticks, ticklabs=ticklabs)
        # keep track of this GeoAxes object for later
        rmd.ax = ax
        return pm

def plot_radar(radarFileName, saveFileName=None, isPreviewRes=False, plotRadius=160, rangeRingStep=None, plot_radial=None, plot_damage=False, plot_lightning=False):
    px = 1/plt.rcParams["figure.dpi"]
    basePath = path.join(getcwd(), "output")
    radarDataDir = path.join(getcwd(), "radarData")
    radarFilePath = path.join(radarDataDir, radarFileName)
    try:
        radar = pyart.io.read(radarFilePath)
    except Exception as e:
        warningString = str(dt.utcnow())+" error reading "+radarFileName+": "+str(e)+"\n"
        logFile = open("warnings.log", "a")
        logFile.write(warningString)
        logFile.close()
        return
    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()})
    if isPreviewRes:
        fig.set_size_inches(768*px, 768*px)
    else:
        fig.set_size_inches(1920*px, 1080*px)
        
    norm, cmap = ctables.registry.get_with_steps("NWSReflectivity", 5, 5)
    cmap.set_under("#00000000")
    cmap.set_over("black")
    ADRADMapDisplay = pyart.graph.RadarMapDisplay(radar)
    plotHandle = plot_ppi_map_modified(ADRADMapDisplay, "reflectivity", 0, resolution="10m", embelish=False, cmap=cmap, norm=norm, colorbar_flag=False, width=2*plotRadius*1000, height=2*plotRadius*1000)
    ADRADMapDisplay.set_aspect_ratio(1)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        ADRADMapDisplay.plot_range_rings(range(0, plotRadius+1, rangeRingStep), col="gray", ls="dotted")
    ax.add_feature(USCOUNTIES.with_scale("5m"), edgecolor="gray")
    if plot_radial is not None:
        ax.plot([radar.longitude["data"][0], radar.longitude["data"][0]+5*np.sin(np.deg2rad(plot_radial))], [radar.latitude["data"][0], radar.latitude["data"][0]+5*np.cos(np.deg2rad(plot_radial))], color="black", linewidth=3)
    radarScanDT = pyart.util.datetime_from_radar(radar)
    radarScanDT = dt(radarScanDT.year, radarScanDT.month, radarScanDT.day, radarScanDT.hour, radarScanDT.minute, radarScanDT.second, tzinfo=timezone.utc)
    if plot_damage:
        damageReports = pd.read_csv("damagereports.csv")
        damageReports["datetimes"] = [dt.strptime(str(damageReports["BEGIN_DATE"][i])+" "+str(damageReports["BEGIN_TIME"][i]), "%m/%d/%y %H%M").replace(tzinfo=timezone.utc) + timedelta(hours=6) for i in range(0, len(damageReports))]
        damageReports = damageReports.set_index(["datetimes"])
        startSearch = radarScanDT - timedelta(hours=1)
        damageReports = damageReports.loc[startSearch:radarScanDT]
        ax.scatter(damageReports["BEGIN_LON"], damageReports["BEGIN_LAT"], s=10*damageReports["MAGNITUDE"], c="black")
    if plot_lightning:
        if path.isdir(path.join(getcwd(), "flashData/2d/")):
            inputDir = path.join(getcwd(), "flashData/2d/")
            for ltgFile in sorted(listdir(inputDir)):
                ltgFile = path.join(inputDir, ltgFile)
                ltgFlash = xr.open_dataset(ltgFile)
                times = ltgFlash["time"].values
                times = [pd.Timestamp(flashTime).to_pydatetime().replace(tzinfo=timezone.utc) for flashTime in times]
                if times[0] > radarScanDT - timedelta(minutes=10):
                    targetIdx = 0
                    for flashTime in times:
                        if flashTime < radarScanDT:
                            targetIdx = times.index(flashTime)
                    lats = ltgFlash["latitude"].values
                    lons = ltgFlash["longitude"].values
                    flashExtent = ltgFlash["flash_extent"].sel(ntimes=targetIdx)
                    ax.contour(lons, lats, flashExtent, cmap="Greys", transform=ccrs.PlateCarree())
                    break
                else:
                    continue
        else:
            inputDir = path.join(getcwd(), "sourceinput")
            for ltgFile in sorted(listdir(inputDir)):
                ltgSrc = readlma.lmafile(path.join(inputDir, ltgFile))
                fileStartTime = ltgSrc.starttime
                fileStartTime = dt(fileStartTime.year, fileStartTime.month, fileStartTime.day, fileStartTime.hour, (fileStartTime.minute - fileStartTime.minute%10), 0, tzinfo=timezone.utc)
                if fileStartTime > radarScanDT - timedelta(minutes=10):
                    with warnings.catch_warnings():
                        warnings.filterwarnings("ignore")
                        ltgData = ltgSrc.readfile()
                    if ltgData.empty:
                        pass
                    else:
                        ltgData["dtobjs"] = [dt(time.year, time.month, time.day, time.hour, time.minute, time.second, tzinfo=timezone.utc) for time in ltgData["Datetime"]]
                        dataToPlot = ltgData.loc[ltgData.dtobjs <= radarScanDT]
                        dataToPlot = dataToPlot.loc[dataToPlot.dtobjs >= radarScanDT - timedelta(minutes=2)]
                        ax.scatter(dataToPlot["lon"], dataToPlot["lat"], s=0.01, c="#00000099", transform=ccrs.PlateCarree())       
    infoString = str()
    if "instrument_name" in radar.metadata.keys():
        insStr = radar.metadata["instrument_name"]
        try:
            insStr = insStr.decode()
        except (UnicodeDecodeError, AttributeError, TypeError):
            pass
        infoString = insStr
    if "sigmet_task_name" in radar.metadata.keys():
        infoString = infoString + " " +radar.metadata["sigmet_task_name"].decode().replace("  ", "")
    elif "vcp_pattern" in radar.metadata.keys():
        infoString = infoString + " VCP-" +str(radar.metadata["vcp_pattern"])
    infoString = infoString + " PPI and Houston LMA VHF **FLASH**\n"
    if "prt" in radar.instrument_parameters:
        prf = np.round(1/np.mean(radar.instrument_parameters["prt"]["data"]), 0)
        infoString = infoString + "Avg. PRF: " + str(prf) + " Hz"
    elevation = np.round(radar.fixed_angle["data"][0], 1)
    infoString = infoString + "    Elevation: " + str(elevation) + "°"
    if "unambiguous_range" in radar.instrument_parameters:
        maxRange = np.round(np.max(radar.instrument_parameters["unambiguous_range"]["data"])/1000, 0)
        infoString = infoString + "    Max Range: " + str(maxRange) + " km\n"
    infoString = infoString + radarScanDT.strftime("%d %b %Y %H:%M:%S UTC")
    ax.set_title(infoString)
    cbax = fig.add_axes([ax.get_position().x0, 0.075, (ax.get_position().width/3), .02])
    fig.colorbar(plotHandle, cax=cbax, orientation="horizontal", extend="neither")
    cbax.set_xlabel("Reflectivity (dBZ)")
    lax = fig.add_axes([ax.get_position().x0+2*(ax.get_position().width/3), 0.015, (ax.get_position().width/3), .1])
    lax.set_aspect(2821/11071)
    plt.setp(lax.spines.values(), visible=False)
    lax.tick_params(left=False, labelleft=False)
    lax.tick_params(bottom=False, labelbottom=False)
    lax.set_xlabel("Plot by Sam Gardner")
    atmoLogo = mpimage.imread("assets/atmoLogo.png")
    lax.imshow(atmoLogo)
    if  isPreviewRes:
        return fig
    else:
        if saveFileName is not None:
            fig.savefig(saveFileName, bbox_inches="tight")
        else:
            fig.savefig(path.join(basePath, str(sorted(listdir(radarDataDir)).index(radarFileName))+".png"), bbox_inches="tight")
        plt.close(fig)

if __name__ == "__main__":
    from itertools import repeat
    radarDataDir = path.join(getcwd(), "radarData")
    with mp.Pool(processes=12) as pool:
        pool.starmap(plot_radar, zip(sorted(listdir(radarDataDir)), repeat(None), repeat(False), repeat(200), repeat(50), repeat(None), repeat(False), repeat(True)))
