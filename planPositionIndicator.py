#!/usr/bin/env python3
# Next-gen HDWX radar plotting script
# Created 7 July 2021 by Sam Gardner <stgardner4@tamu.edu>

from datetime import datetime as dt
from datetime import timedelta, timezone
from netCDF4 import Dataset
import pyart
from matplotlib import pyplot as plt
from os import path, listdir
from cartopy import crs as ccrs
from metpy.plots import ctables
from metpy.plots import USCOUNTIES
import numpy as np
import warnings
import multiprocessing as mp
from matplotlib import image as mpimage
import pandas as pd
from pyxlma.lmalib.io import read as lma_read
from pyxlma.lmalib.flash.cluster import cluster_flashes
from pyxlma.lmalib.grid import  create_regular_grid, assign_regular_bins, events_to_grid

basePath = path.dirname(path.realpath(__file__))
axExtent = [-98.5, -95, 30.25, 32.25]

def plot_radar(radarFileName, saveFileName=None, isPreviewRes=False, plotRadius=160, rangeRingStep=None, plot_radial=None, plot_damage=False, plot_lightning=False, field="reflectivity"):
    px = 1/plt.rcParams["figure.dpi"]
    outputPath = path.join(basePath, "output")
    radarDataDir = path.join(basePath, "radarData")
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
    ADRADMapDisplay = pyart.graph.RadarMapDisplay(radar)
    if field == "reflectivity":
        norm, cmap = ctables.registry.get_with_steps("NWSReflectivity", 5, 5)
        cmap.set_under("#00000000")
        cmap.set_over("black")
        cbStr = "Reflectivity (dBZ)"
        ADRADMapDisplay.plot_ppi_map("reflectivity", 0, resolution="10m", embelish=False, cmap=cmap, norm=norm, colorbar_flag=False, width=2*plotRadius*1000, height=2*plotRadius*1000)
    elif field == "velocity":
        norm, cmap = ctables.registry.get_with_steps("NWS8bitVel", -100, 1)
        cbStr = "Velocity (m/s)"
        ADRADMapDisplay.plot_ppi_map("velocity", 1, resolution="10m", embelish=False, cmap=cmap, norm=norm, colorbar_flag=False, width=2*plotRadius*1000, height=2*plotRadius*1000)
    elif field == "differential_reflectivity":
        cbStr = "Differential Reflectivity (dB)"
        ADRADMapDisplay.plot_ppi_map("differential_reflectivity", 0, resolution="10m", embelish=False, cmap="pyart_RefDiff", vmin=-1, vmax=8, colorbar_flag=False, width=2*plotRadius*1000, height=2*plotRadius*1000)
    plotHandle = ax.get_children()[0]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        ADRADMapDisplay.plot_range_rings(range(0, plotRadius+1, rangeRingStep), col="gray", ls="dotted")
    ax.add_feature(USCOUNTIES.with_scale("5m"), edgecolor="gray", zorder=3)
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
    ltgType = False
    flashContours = ""
    if plot_lightning:
        ltgFiles = sorted(listdir(path.join(basePath, "lightningin")))
        radarFiles = sorted(listdir(path.join(basePath, "radarData")))
        i_want_exact = False
        if "0600.dat.gz" in ltgFiles[0]:
            numMins = 10
            targetTimes = [dt(radarScanDT.year, radarScanDT.month, radarScanDT.day, radarScanDT.hour, (radarScanDT.minute) - radarScanDT.minute % 10, 0).replace(tzinfo=None), radarScanDT.replace(tzinfo=None)]
        elif "0060.dat.gz" in ltgFiles[0]:
            numMins = 1
            lastRadFileName = radarFiles[radarFiles.index(radarFileName) - 1]
            if "TAMU" in lastRadFileName:
                radFileArr = lastRadFileName.split("_")
                timeOfLastScan = dt.strptime(radFileArr[1]+radFileArr[2], "%Y%m%d%H%M")
            elif "V06" in lastRadFileName:
                radFileArr = lastRadFileName[4:].split("_")
                timeOfLastScan = dt.strptime(radFileArr[0]+radFileArr[1], "%Y%m%d%H%M%S")
            targetTimes = [timeOfLastScan.replace(tzinfo=None), radarScanDT.replace(tzinfo=None)]
        targetLtgFiles = list()
        for file in ltgFiles:
            timeOfFileArr = file.split("_")
            ltgFileTime = dt.strptime("20"+timeOfFileArr[1]+timeOfFileArr[2], "%Y%m%d%H%M%S").replace(tzinfo=None) + timedelta(minutes=1)
            if i_want_exact:
                if dt(ltgFileTime.year, ltgFileTime.month, ltgFileTime.day, ltgFileTime.hour, ltgFileTime.minute, 0) == dt(radarScanDT.year, radarScanDT.month, radarScanDT.day, radarScanDT.hour, radarScanDT.minute, 0):
                    targetLtgFiles = [path.join(basePath, "lightningin", file)]
            else:
                if ltgFileTime >= targetTimes[0] and ltgFileTime < targetTimes[1]:
                    targetLtgFiles.append(path.join(basePath, "lightningin", file))
        if len(targetLtgFiles) > 0:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                # Read in LMA data
                lmaData, startTimeOfPlot = lma_read.dataset(targetLtgFiles)
            timeOfPlot = startTimeOfPlot + numMins*timedelta(minutes=len(targetLtgFiles))
            dttuple = (np.datetime64(startTimeOfPlot), np.datetime64(timeOfPlot))
            grid_dt = np.asarray(60, dtype='m8[s]')
            grid_t0 = np.asarray(dttuple[0]).astype('datetime64[ns]')
            grid_t1 = np.asarray(dttuple[1]).astype('datetime64[ns]')
            time_range = (grid_t0, grid_t1+grid_dt, grid_dt)
            # We only want events with chi^2 less than 1
            lmaData = lmaData[{"number_of_events":(lmaData.event_chi2 <= 1.0)}]
            lmaData = cluster_flashes(lmaData)
            lat_range = (axExtent[2], axExtent[3], 0.025)
            lon_range = (axExtent[0], axExtent[1], 0.025)
            alt_range = (0, 18e3, 1.0e3)
            grid_edge_ranges ={
                'grid_latitude_edge':lat_range,
                'grid_longitude_edge':lon_range,
                'grid_altitude_edge':alt_range,
                'grid_time_edge':time_range,
            }
            grid_center_names ={
                'grid_latitude_edge':'grid_latitude',
                'grid_longitude_edge':'grid_longitude',
                'grid_altitude_edge':'grid_altitude',
                'grid_time_edge':'grid_time',
            }
            event_coord_names = {
                'event_latitude':'grid_latitude_edge',
                'event_longitude':'grid_longitude_edge',
                'event_altitude':'grid_altitude_edge',
                'event_time':'grid_time_edge',
            }
            grid_ds = create_regular_grid(grid_edge_ranges, grid_center_names)
            ds_ev = assign_regular_bins(grid_ds, lmaData, event_coord_names, pixel_id_var="event_pixel_id", append_indices=True)
            grid_spatial_coords=['grid_time', None, 'grid_latitude', 'grid_longitude']
            event_spatial_vars = ('event_altitude', 'event_latitude', 'event_longitude')
            griddedLmaData = events_to_grid(ds_ev, grid_ds, min_points_per_flash=3, pixel_id_var="event_pixel_id", event_spatial_vars=event_spatial_vars, grid_spatial_coords=grid_spatial_coords)
            griddedLmaData = griddedLmaData.isel(grid_time=0)
            try:
                flashContours = ax.contourf(griddedLmaData.flash_extent_density.grid_longitude, griddedLmaData.flash_extent_density.grid_latitude, griddedLmaData.flash_extent_density.data, levels=np.arange(1, 14.01, 0.1), cmap="plasma", alpha=0.5, transform=ccrs.PlateCarree())
            except Exception as e:
                if "GEOSContains" in str(e):
                    return
                else:
                    raise e
            
    infoString = str()
    if path.exists("stormlocations.csv"):
        stormLocs = pd.read_csv("stormlocations.csv")
        infoString = "Storm-Following "
        stormLocs["londiff"] = stormLocs["lonmax"] - stormLocs["lonmin"]
        maxLonDiff = max(stormLocs["londiff"])
        stormLocs["latdiff"] = stormLocs["latmax"] - stormLocs["latmin"]
        maxLatDiff = max(stormLocs["latdiff"])
        stormLoc = stormLocs.loc[stormLocs["time"] == radarScanDT.strftime("%d.%H.%M.%S")]
        minLon = np.mean([float(stormLoc["lonmin"]), float(stormLoc["lonmax"])]) - maxLonDiff/2
        maxLon = np.mean([float(stormLoc["lonmin"]), float(stormLoc["lonmax"])]) + maxLonDiff/2
        minLat = np.mean([float(stormLoc["latmin"]), float(stormLoc["latmax"])]) - maxLatDiff/2
        maxLat = np.mean([float(stormLoc["latmin"]), float(stormLoc["latmax"])]) + maxLatDiff/2
        try:
            ax.set_extent([minLon, maxLon, minLat, maxLat])
        except Exception as e:
            print("ERR! ERR! ERR!")
            print(radarScanDT)
    else:
        ax.set_extent(axExtent)
        infoString = ""
    if "instrument_name" in radar.metadata.keys():
        insStr = radar.metadata["instrument_name"]
        try:
            insStr = insStr.decode()
        except (UnicodeDecodeError, AttributeError, TypeError):
            pass
        infoString = infoString+insStr
    if "sigmet_task_name" in radar.metadata.keys():
        infoString = infoString + " " +radar.metadata["sigmet_task_name"].decode().replace("  ", "")
    elif "vcp_pattern" in radar.metadata.keys():
        infoString = infoString + " VCP-" +str(radar.metadata["vcp_pattern"])
    if type(ltgType) == type("THIS IS A STRING"):
        infoString = infoString+" PPI and Houston LMA "+ltgType+"\n"
    else:
        infoString = infoString+" PPI\n"
    if "prt" in radar.instrument_parameters:
        prf = np.round(1/np.mean(radar.instrument_parameters["prt"]["data"]), 0)
        infoString = infoString + "Avg. PRF: " + str(prf) + " Hz"
    elevation = np.round(radar.fixed_angle["data"][0], 1)
    infoString = infoString + "    Elevation: " + str(elevation) + "°"
    if "unambiguous_range" in radar.instrument_parameters:
        maxRange = np.round(np.max(radar.instrument_parameters["unambiguous_range"]["data"])/1000, 0)
        infoString = infoString + "    Max Range: " + str(maxRange) + " km\n"
    infoString = infoString + radarScanDT.strftime("%d %b %Y %H:%M:%S UTC")
    ax.gridlines(draw_labels=True)
    ax.set_position([0.05, 0.11, .9, .87])
    cbax = fig.add_axes([0, 0, (ax.get_position().width/3), .025])
    fig.colorbar(plotHandle, cax=cbax, orientation="horizontal", extend="neither")
    cbax.set_position([0.05, ax.get_position().y0-.1-cbax.get_position().height, cbax.get_position().width, cbax.get_position().height])
    if flashContours is not "":
        cbax2 = fig.add_axes([0, 0, (ax.get_position().width/3), .025])
        fig.colorbar(flashContours, cax=cbax2, orientation="horizontal", label="Flash Extent Density", extend="max").set_ticks(np.arange(1, 14.01, 1))
        cbax2.set_position([0.05, ax.get_position().y0-.1-cbax.get_position().height-.01-cbax2.get_position().height, cbax2.get_position().width, cbax2.get_position().height])
    tax = fig.add_axes([0,0,(ax.get_position().width/3),.05])
    tax.text(0.5, 0.5, infoString, horizontalalignment="center", verticalalignment="center", fontsize=16)
    tax.set_xlabel("Plot by Sam Gardner")
    plt.setp(tax.spines.values(), visible=False)
    tax.tick_params(left=False, labelleft=False)
    tax.tick_params(bottom=False, labelbottom=False)
    tax.set_position([(ax.get_position().width/2)-(tax.get_position().width/2)+ax.get_position().x0,ax.get_position().y0-.1-tax.get_position().height,tax.get_position().width,tax.get_position().height], which="both")
    lax = fig.add_axes([0,0,(ax.get_position().width/3),1])
    lax.set_aspect(2821/11071)
    lax.axis("off")
    plt.setp(lax.spines.values(), visible=False)
    atmoLogo = mpimage.imread("assets/atmoLogo.png")
    lax.imshow(atmoLogo)
    lax.set_position([.95-lax.get_position().width, ax.get_position().y0-.1-lax.get_position().height, lax.get_position().width, lax.get_position().height], which="both")
    fig.set_facecolor("white")
    if  isPreviewRes:
        fig.set_size_inches(768*px, 768*px)
        return fig
    else:
        fig.set_size_inches(2560*px, 1440*px)
        if saveFileName is not None:
            fig.savefig(saveFileName)
        else:
            fig.savefig(path.join(outputPath, str(sorted(listdir(radarDataDir)).index(radarFileName))+".png"))
        plt.close(fig)

if __name__ == "__main__":
    from itertools import repeat
    radarDataDir = path.join(basePath, "radarData")
    # plot_radar(sorted(listdir(radarDataDir))[0], None, False, 200, 50, None, False, True, "reflectivity")
    # exit()
    with mp.Pool(processes=8) as pool:
        pool.starmap(plot_radar, zip(sorted(listdir(radarDataDir)), repeat(None), repeat(False), repeat(200), repeat(50), repeat(None), repeat(False), repeat(True), repeat("reflectivity")))
