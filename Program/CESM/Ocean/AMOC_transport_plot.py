#Program plots the AMOC strength

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
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker


#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	transport	= fh.variables['Transport'][:]	#AMOC strength (Sv)

	fh.close()

	return time, transport

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

depth_min 	= 0
depth_max	= 1000

time, transport			    = ReadinData(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-100, 2500], 16, 19, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot([0.1/0.0003, 0.1/0.0003], [-5, 20], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.2/0.0003, 0.2/0.0003], [-5, 20], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.3/0.0003, 0.3/0.0003], [-5, 20], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.4/0.0003, 0.4/0.0003], [-5, 20], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.5/0.0003, 0.5/0.0003], [-5, 20], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.6/0.0003, 0.6/0.0003], [-5, 20], linestyle = '--', color = 'c', linewidth = 1)

ax.text(0.1/0.0003, 20, '0.1 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.2/0.0003, 20, '0.2 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.3/0.0003, 20, '0.3 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.4/0.0003, 20, '0.4 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.5/0.0003, 20, '0.5 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.6/0.0003, 20, '0.6 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)

graph_control	= plot(time, transport, '-k', linewidth = 0.5)

ax.set_xlabel('Model year')
ax.set_ylabel('Volume transport (Sv)')
ax.set_xlim(1, 2200)
ax.set_ylim(-2, 22)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])
ax.set_title('a) AMOC strength at 26$^{\circ}$N')

ax.quiver(1758, 14, 0, -0.4, scale = 5, color = 'r', zorder = 10, width = 0.005)
ax.text(1758, 14.05, 'AMOC tipping\n(1758)', verticalalignment='bottom', horizontalalignment='center', color = 'r', fontsize = 10)

ax.plot([1, 50], [12.7, 12.7], '-b')
ax.plot([1, 1], [12.3, 13.1], '-b', zorder = 100, clip_on = False)
ax.plot([50, 50], [12.3, 13.1], '-b')
ax.text(30, 11.05, '1 - 50', verticalalignment='center', horizontalalignment='center', rotation = 'vertical', color = 'b', fontsize = 8)

ax.plot([1701, 1750], [6, 6], '-b')
ax.plot([1701, 1701], [5.6, 6.4], '-b')
ax.plot([1750, 1750], [5.6, 6.4], '-b')
ax.text(1730, 3.2, '1701 - 1750', verticalalignment='center', horizontalalignment='center', rotation = 'vertical', color = 'b', fontsize = 8)

ax.plot([2151, 2200], [3.7, 3.7], '-b')
ax.plot([2151, 2151], [3.3, 4.1], '-b')
ax.plot([2200, 2200], [3.3, 4.1], '-b', zorder = 100, clip_on = False)
ax.text(2178, 6.5, '2151 - 2200', verticalalignment='center', horizontalalignment='center', rotation = 'vertical', color = 'b', fontsize = 8)

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.06, 0.12, 0.45, 0.40], projection = ccrs.Orthographic(-30, 10))

ax2.coastlines(resolution='50m')
ax2.gridlines()
ax2.add_feature(cfeature.LAND, zorder=10)
ax2.set_global()


lon     = np.arange(-1, 360)
lat     = np.arange(-90, 91)
field   = np.ones((len(lat), len(lon))) * -0.35
CS      = ax2.contourf(lon, lat, field, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())


lon     = np.arange(-100, -5)
lat     = np.arange(20, 43)
field   = np.ones((len(lat), len(lon))) * 0.35
CS      = ax2.contourf(lon, lat, field, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

lon     = np.arange(-100, 3)
lat     = np.arange(42, 51)
field   = np.ones((len(lat), len(lon))) * 0.35
CS      = ax2.contourf(lon, lat, field, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

ax2.text(320, 38, '$+F_H$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=11, transform=ccrs.PlateCarree())
ax2.text(340, -10, '$-F_H$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=11, transform=ccrs.PlateCarree())

x_1	= np.arange(-81, -9.99, 0.1)
y_1	= np.zeros(len(x_1)) + 26.0
y_2	= np.arange(24, 28.01, 0.1)
x_2	= np.zeros(len(y_2)) + x_1[0]
y_3	= np.arange(24, 28.01, 0.1)
x_3	= np.zeros(len(y_3)) + x_1[-1]

ax2.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_2, y_2, '-k', linewidth = 2.0
, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

x_1	= np.arange(-60, 20.01, 0.1)
y_1	= np.zeros(len(x_1)) - 34
y_2	= np.arange(-37, -30.99, 0.1)
x_2	= np.zeros(len(y_2)) + x_1[0]
y_3	= np.arange(-37, -30.99, 0.1)
x_3	= np.zeros(len(y_3)) + x_1[-1]

ax2.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_2, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

#-----------------------------------------------------------------------------------------

show()