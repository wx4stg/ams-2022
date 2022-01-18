#!/usr/bin/env python3
# Next-gen HDWX radar cross-section visualization
# Created 18 August 2021 by Sam Gardner <stgardner4@tamu.edu>

from datetime import datetime as dt, timedelta, timezone
from turtle import distance
from matplotlib import pyplot as plt
from metpy.plots import ctables
from cartopy import geodesic
from os import path, getcwd, listdir
import pyart
import numpy as np
from matplotlib import image as mpimage
import warnings
from pyxlma.lmalib.io import read as lma_read
from pyxlma.lmalib.flash.cluster import cluster_flashes
from pyxlma.lmalib.grid import  create_regular_grid, assign_regular_bins, events_to_grid

def plot_crosssection(radarFileName, saveFileName=None, requestedLatitude=0, isPreviewRes=False, plotRadius=160, rangeRingStep=10, plot_lightning=True):
    px = 1/plt.rcParams["figure.dpi"]
    basePath = path.dirname(path.realpath(__file__))
    outDir = path.join(getcwd(), "output")
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
    fig, ax = plt.subplots(1, 1)
    if isPreviewRes:
        fig.set_size_inches(768*px, 768*px)
    else:
        fig.set_size_inches(1920*px, 1080*px)
    norm, cmap = ctables.registry.get_with_steps("NWSReflectivity", 5, 5)
    cmap.set_under("#00000000")
    cmap.set_over("black")
    ADRADGrid = pyart.map.grid_from_radars(radar, grid_shape=(21, 321, 321), grid_limits=((0, 20000), (-160000, 160000), (-160000, 160000)))
    ADRADDisplay = pyart.graph.GridMapDisplay(ADRADGrid)
    ADRADDisplay.plot_latitude_slice("reflectivity", lat=requestedLatitude, cmap=cmap, norm=norm, colorbar_flag=False, title_flag=False, ax=ax)
    plotHandle = ax.get_children()[0]
    radarScanDT = pyart.util.datetime_from_radar(radar)
    radarScanDT = dt(radarScanDT.year, radarScanDT.month, radarScanDT.day, radarScanDT.hour, radarScanDT.minute, radarScanDT.second, tzinfo=timezone.utc)
    flashContours = ""
    if plot_lightning:
        axExtent = [-98.5, -95, 30.25, 32.25]
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
            lat_range = (axExtent[2], axExtent[3], 0.0125)
            lon_range = (axExtent[0], axExtent[1], 0.0125)
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
            grid_spatial_coords=['grid_time', 'grid_altitude', 'grid_latitude', 'grid_longitude']
            event_spatial_vars = ('event_altitude', 'event_latitude', 'event_longitude')
            griddedLmaData = events_to_grid(ds_ev, grid_ds, min_points_per_flash=3, pixel_id_var="event_pixel_id", event_spatial_vars=event_spatial_vars, grid_spatial_coords=grid_spatial_coords)
            griddedLmaData = griddedLmaData.flash_extent_density.sum("grid_time")
            griddedLmaData = griddedLmaData.sel(grid_latitude=slice(requestedLatitude-.03, requestedLatitude+.03))
            lmaLongs = griddedLmaData.grid_longitude.data
            radarLong = radar.longitude["data"][0]
            distanceFromRadar = np.array([geodesic.Geodesic().inverse([radarLong, requestedLatitude], [gridLong, requestedLatitude])[0][0]/1000 for gridLong in lmaLongs])
            distanceFromRadar = np.tile(np.array([distanceFromRadar.transpose()]), [griddedLmaData.grid_altitude.data.shape[0], 1])
            heightAboveGround = np.tile(griddedLmaData.grid_altitude.data/1000, [distanceFromRadar.shape[1], 1]).transpose()
            griddedLmaData = griddedLmaData.sum("grid_latitude")
            zeroMask = np.where(griddedLmaData.data < 1.0, 1, 0)
            print(np.amax(griddedLmaData.data))
            try:
                flashContours = ax.pcolormesh(distanceFromRadar, heightAboveGround, np.ma.masked_array(griddedLmaData.data, mask=zeroMask), vmin=1, vmax=75, cmap="plasma", alpha=0.75)
            except Exception as e:
                print(e)
                if "GEOSContains" in str(e):
                    return
                else:
                    raise e

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
    infoString = infoString + " Cross-Section\n"
    if "prt" in radar.instrument_parameters:
        prf = np.round(1/np.mean(radar.instrument_parameters["prt"]["data"]), 0)
        infoString = infoString + "Avg. PRF: " + str(prf) + " Hz"
    elevation = np.round(radar.fixed_angle["data"][0], 1)
    infoString = infoString + "    Elevation: " + str(elevation) + "°"
    if "unambiguous_range" in radar.instrument_parameters:
        maxRange = np.round(np.max(radar.instrument_parameters["unambiguous_range"]["data"])/1000, 0)
        infoString = infoString + "    Max Range: " + str(maxRange) + " km\n"
    infoString = infoString + pyart.util.datetime_from_radar(radar).strftime("%d %b %Y %H:%M:%S UTC")
    ax.set_xlim([0, plotRadius])
    ax.set_xticks(range(0, plotRadius+1, rangeRingStep))
    ax.set_ylim([0, 20])
    ax.set_title(infoString)
    cbax = fig.add_axes([ax.get_position().x0, 0.05, (ax.get_position().width/3), .02])
    fig.colorbar(plotHandle, cax=cbax, orientation="horizontal", extend="neither")
    cbax.set_xlabel("Reflectivity (dBZ)")
    if flashContours != "":
        cbax2 = fig.add_axes([0, 0, (ax.get_position().width/3), .025])
        fig.colorbar(flashContours, cax=cbax2, orientation="horizontal", label="Flash Extent Density", extend="neither").set_ticks([1]+list(np.arange(5, 75.1, 5)))
        cbax2.set_position([cbax.get_position().x0, ax.get_position().y0-.1-cbax.get_position().height-.01-cbax2.get_position().height, cbax2.get_position().width, cbax2.get_position().height])
    lax = fig.add_axes([ax.get_position().x0+2*(ax.get_position().width/3), 0, (ax.get_position().width/3), .1])
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
            fig.savefig("rhitest.png", bbox_inches="tight")
        plt.close(fig)
