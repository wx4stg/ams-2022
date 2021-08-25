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


def plot_azimuth_to_rhi_modified(
        rd, field, target_azimuth, mask_tuple=None,
        vmin=None, vmax=None, norm=None, cmap=None, mask_outside=False,
        title=None, title_flag=True,
        axislabels=(None, None), axislabels_flag=True,
        colorbar_flag=True, colorbar_label=None,
        colorbar_orient='vertical', edges=True, gatefilter=None,
        reverse_xaxis=None, filter_transitions=True,
        ax=None, fig=None, ticks=None, ticklabs=None,
        raster=False, **kwargs):
    # parse parameters
    ax, fig = pyart.graph.common.parse_ax_fig(ax, fig)
    vmin, vmax = pyart.graph.common.parse_vmin_vmax(rd._radar, field, vmin, vmax)
    cmap = pyart.graph.common.parse_cmap(cmap, field)

    data, x, y, z = rd._get_azimuth_rhi_data_x_y_z(
        field, target_azimuth, edges, mask_tuple,
        filter_transitions, gatefilter)

    # mask the data where outside the limits
    if mask_outside:
        data = np.ma.masked_invalid(data)
        data = np.ma.masked_outside(data, vmin, vmax)

    # plot the data
    R = np.sqrt(x ** 2 + y ** 2) * np.sign(y)
    if reverse_xaxis is None:
        # reverse if all distances (nearly, up to 1 m) negative.
        reverse_xaxis = np.all(R < 1.)
    if reverse_xaxis:
        R = -R
    if norm is not None:  # if norm is set do not override with vmin/vmax
        vmin = vmax = None
    pm = ax.pcolormesh(
        R, z, data, vmin=vmin, vmax=vmax, cmap=cmap, norm=norm, **kwargs)

    if raster:
        pm.set_rasterized(True)

    if title_flag:
        rd._set_az_rhi_title(field, target_azimuth, title, ax)

    if axislabels_flag:
        rd._label_axes_rhi(axislabels, ax)

    # add plot and field to lists
    rd.plots.append(pm)
    rd.plot_vars.append(field)

    if colorbar_flag:
        rd.plot_colorbar(
            mappable=pm, label=colorbar_label, orient=colorbar_orient,
            field=field, ax=ax, fig=fig, ticks=ticks, ticklabs=ticklabs)
    return pm

def plot_crosssection(radarFileName, saveFileName=None, azimuth=0, isPreviewRes=False, plotRadius=160, rangeRingStep=10):
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
    ADRADDisplay = pyart.graph.RadarDisplay(radar)
    plotHandle = plot_azimuth_to_rhi_modified(ADRADDisplay, "reflectivity", azimuth, norm=norm, cmap=cmap, colorbar_flag=False)
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
    infoString = infoString + " PPI\n"
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
