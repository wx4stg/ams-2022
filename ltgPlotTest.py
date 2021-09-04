import xarray as xr
from matplotlib import pyplot as plt
from cartopy import crs as ccrs



fig, ax = plt.subplots(1, 1, subplot_kw={'projection': ccrs.PlateCarree()})
ltgFlash = xr.open_dataset("flashData/2d/HSTN_20210409_021000_600_10src_0.0103deg-dx_flash_extent.nc")
lon = ltgFlash.variables["longitude"].values
lat = ltgFlash.variables["latitude"].values
flashEx = ltgFlash["flash_extent"].values
ax.contour(lon, lat, flashEx[0,:,:])
ax.gridlines(draw_labels=True)
fig.savefig("test.png")