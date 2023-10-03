#Program plots the Southern hemispheric sea-ice distribution

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

month_ice	= 9
#-----------------------------------------------------------------------------------------

fh 	= netcdf.Dataset(directory+'Data/Antarctic_sea_ice_fraction_September/CESM_year_0001-0050.nc', 'r')

lon		= fh.variables['lon'][:]
lat		= fh.variables['lat'][:] 			
fraction_1	= np.mean(fh.variables['Fraction'][:], axis = 0)

fh.close()

fh 	= netcdf.Dataset(directory+'Data/Antarctic_sea_ice_fraction_September/CESM_year_2151-2200.nc', 'r')
			
fraction_2	= np.mean(fh.variables['Fraction'][:], axis = 0)

fh.close()

#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Orthographic(0, -90)})

CS      = ax.contourf(lon, lat, fraction_2, levels = np.arange(15, 100.01, 4.25), cmap = 'tab20b', transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon, lat, fraction_1, levels = [15], colors = 'darkblue', linewidths = 2, transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = [15, 40, 60, 80, 100], cax=ax_cb)
cbar.set_label('Sea-ice fraction ($\%$)')

ax.gridlines(zorder=10)
ax.add_feature(cfeature.LAND, zorder=10)
ax.coastlines()

ax.scatter([-90, 0, 90, 180], [-35, -35, -35, -35], marker = 'o', edgecolor = 'w' , s = 6, facecolors='none', transform=ccrs.PlateCarree())

graphs	= ax.plot([50, 50], [-90, -90], '-', color = 'darkblue', linewidth = 2.0, label = '1 - 50', transform=ccrs.PlateCarree(), zorder =0)

legend_labels = [l.get_label() for l in graphs]

ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1).set_zorder(12)


ax.set_title('b) Antarctic sea-ice extent, September (2151 - 2200)')

show()



