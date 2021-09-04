from matplotlib import pyplot as plt
from cartopy import crs as ccrs
from datetime import datetime as dt
from datetime import timedelta
from os import path, listdir
from netCDF4 import Dataset
import numpy as np


def centers_to_edges(x):
    xedge=np.zeros(x.shape[0]+1)
    xedge[1:-1] = (x[:-1] + x[1:])/2.0
    dx = np.mean(np.abs(xedge[2:-1] - xedge[1:-2]))
    xedge[0] = xedge[1] - dx
    xedge[-1] = xedge[-2] + dx
    return xedge

if __name__ == "__main__":
    radarScanDt = dt(2021, 4, 9, 1, 52, 38, 0)
    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()})
    inputDir = "flashData/2d"
    for ltgFile in sorted(listdir(inputDir)):
        ltgFile = path.join(inputDir, ltgFile)
        ltgFlash = Dataset(ltgFile)
        flashData = ltgFlash.variables
        times = flashData["time"]
        base_date = dt.strptime(times.units, "seconds since %Y-%m-%d %H:%M:%S")
        time_delta = timedelta(0, float(times[0]),0)
        start_time = base_date + time_delta
        if radarScanDt - timedelta(minutes=10) < start_time:
            flashDims = ltgFlash.dimensions
            lon = flashData["longitude"]
            lat = flashData["latitude"]
            grid = flashData["flash_extent"]
            grid_dims = grid.dimensions
            name_to_idx = dict((k, i) for i, k in enumerate(grid_dims))
            grid_t_idx = name_to_idx[times.dimensions[0]]
            n_frames = times.shape[0]
            density_maxes = []
            total_counts = []
            all_t = []
            xedge = centers_to_edges(lon)
            x_range = xedge.max() - xedge.min()
            yedge = centers_to_edges(lat)
            y_range = yedge.max() - yedge.min()
            dx = (xedge[1]-xedge[0])
            min_count, max_count = 1, grid[:].max()
            if (max_count == 0) | (max_count == 1 ):
                max_count = min_count+1
            default_vmin = -1.0
            if np.log10(max_count) <= default_vmin:
                vmin_count = np.log10(max_count) + default_vmin
            else:
                vmin_count = default_vmin
            indexer = [slice(None),]*len(grid.shape)

            frame_start_times = []
            for i in range(n_frames):
                frame_start = base_date + timedelta(seconds=float(times[i]))
                frame_start_times.append(frame_start)
                indexer[grid_t_idx] = i
                print(frame_start)
                if frame_start > radarScanDt:
                    break
            print(frame_start)
            density = grid[indexer]
            ax.pcolormesh(xedge,yedge, np.log10(density.transpose()), vmin=vmin_count,vmax=np.log10(max_count), cmap="Greys")
            ax.gridlines(draw_labels=True)
            fig.savefig("test"+str(i)+".png")
            break