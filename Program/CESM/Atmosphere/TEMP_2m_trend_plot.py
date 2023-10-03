#Program plots the 2-meter surface temperature trend

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
from scipy import stats
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

fh 		           = netcdf.Dataset(directory+'Atmosphere/TEMP_2m_trend_year_1750-1850.nc', 'r')

lon			        = fh.variables['lon'][:]
lat			        = fh.variables['lat'][:] 			
temp_trend_1		= fh.variables['TEMP_trend'][:] 
temp_trend_1_sig	= fh.variables['TEMP_trend_sig'][:]

fh.close()

#-----------------------------------------------------------------------------------------

fh 			            = netcdf.Dataset(directory+'Atmosphere/TEMP_2m_monthly_trend_year_1750-1850.nc', 'r')

#Select February		
temp_trend_1_month	    = fh.variables['TEMP_trend'][1] 
temp_trend_1_month_sig	= fh.variables['TEMP_trend_sig'][1] 

fh.close()

#-----------------------------------------------------------------------------------------
#Rescale the temperature plot
scale	= 5.0
cut_off	= 5

temp_trend_1[temp_trend_1 < -cut_off]			= (temp_trend_1[temp_trend_1 < -cut_off] - -cut_off) / scale - cut_off
temp_trend_1[temp_trend_1 > cut_off]			= (temp_trend_1[temp_trend_1 > cut_off] - cut_off) / scale + cut_off
temp_trend_1_month[temp_trend_1_month < -cut_off]	= (temp_trend_1_month[temp_trend_1_month < -cut_off] - -cut_off) / scale - cut_off
temp_trend_1_month[temp_trend_1_month > cut_off]	= (temp_trend_1_month[temp_trend_1_month > cut_off] - cut_off) / scale + cut_off

#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Robinson()})

CS      = ax.contourf(lon, lat, temp_trend_1, levels = np.arange(-8, 8.01, 0.5), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = [-8, -6, -4, -2, 0, 2, 4, 6, 8], cax=ax_cb)
cbar.ax.set_yticklabels([-20, -10, -4, -2, 0, 2, 4, 10, 20])
cbar.set_label('2-meter temperature trend ($^{\circ}$C per century)')

ax.set_global()
ax.gridlines()
ax.coastlines()

for lat_i in range(0, len(lat), 3):
	for lon_i in range(0, len(lon), 3):
		#Determine significant difference

		if temp_trend_1_sig[lat_i, lon_i] <= 0.95:
			#Non-significant difference
			ax.scatter(lon[lon_i], lat[lat_i], marker = 'o', edgecolor = 'k' , s = 6, facecolors='none', transform=ccrs.PlateCarree())

ax.set_title('a) Yearly 2-meter temperature trend (1750 - 1850)')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax        = plt.subplots(subplot_kw={'projection': ccrs.Stereographic(central_longitude=0, central_latitude=40)})

CS      = ax.contourf(lon, lat, temp_trend_1_month, levels = np.arange(-12, 12.01, 0.5), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = [-12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12], cax=ax_cb)
cbar.ax.set_yticklabels([-40, -30, -20, -10, -4, -2, 0, 2, 4, 10, 20, 30, 40])
cbar.set_label('2-meter temperature trend ($^{\circ}$C per century)')

ax.coastlines('50m')
gl1	     = ax.gridlines(draw_labels=True, dms = True, x_inline=False, y_inline=False, linewidth = 0.0)
gl1.top_labels = False
gl1.right_labels = False
gl1.xlocator = mticker.FixedLocator([-30, 0, 30])
gl1.ylocator = mticker.FixedLocator([20, 30, 40, 50])
gl1.xlabel_style = {'rotation':0}
gl2 	= ax.gridlines(draw_labels=False, dms = True, x_inline=False, y_inline=False)

ax.set_extent([-41, 41, 25, 78], ccrs.PlateCarree())

for lat_i in range(0, len(lat)):
	for lon_i in range(0, len(lon)):
		#Determine significant difference

		if temp_trend_1_month_sig[lat_i, lon_i] <= 0.95:
			#Non-significant difference
			ax.scatter(lon[lon_i], lat[lat_i], marker = 'o', edgecolor = 'k' , s = 6, facecolors='none', transform=ccrs.PlateCarree())


ax.scatter(-0.1, 51.5, marker = 'o', color = 'r' , s = 12, transform=ccrs.PlateCarree(), zorder = 10)
ax.scatter(-3.7, 40.4, marker = 'o', color = 'r' , s = 12, transform=ccrs.PlateCarree(), zorder = 10)
ax.scatter(16.4, 48.2, marker = 'o', color = 'r' , s = 12, transform=ccrs.PlateCarree(), zorder = 10)
ax.scatter(5.3, 60.4, marker = 'o', color = 'r' , s = 12, transform=ccrs.PlateCarree(), zorder = 10)
ax.scatter(-21.9, 64.1, marker = 'o', color = 'r' , s = 12, transform=ccrs.PlateCarree(), zorder = 10)

ax.set_title('b) February 2-meter temperature trend (1750 - 1850)')

show()
