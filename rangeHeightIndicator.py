#!/usr/bin/env python3
# Next-gen HDWX radar cross-section visualization
# Created 18 August 2021 by Sam Gardner <stgardner4@tamu.edu>

from datetime import datetime as dt
from matplotlib import pyplot as plt
from metpy.plots import ctables
from os import path, getcwd
import pyart
import numpy as np
from matplotlib import image as mpimage


def plot_crosssection(radarFileName, saveFileName=None, requestedLatitude=0, isPreviewRes=False, plotRadius=160, rangeRingStep=10, plotLightning):
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
    ADRADDisplay.plot_latitude_slice("reflectivity", lat=requestedLatitude, cmap=cmap, colorbar_flag=False, title_flag=False, ax=ax)
    plotHandle = ax.get_children()[0]
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
    print()
    ax.set_xlim([0, plotRadius])
    ax.set_xticks(range(0, plotRadius+1, rangeRingStep))
    ax.set_ylim([0, (plotRadius*np.tan(np.deg2rad(max(radar.elevation["data"]))))])
    ax.set_title(infoString)
    cbax = fig.add_axes([ax.get_position().x0, 0.05, (ax.get_position().width/3), .02])
    fig.colorbar(plotHandle, cax=cbax, orientation="horizontal", extend="neither")
    cbax.set_xlabel("Reflectivity (dBZ)")
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
