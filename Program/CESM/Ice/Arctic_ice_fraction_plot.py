#Program plots the Northern hemispheric sea-ice distribution

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

def ConverterField(index_break, field):
	"""Shifts field, where it starts at 0E and ends at 360E"""

	new_field	= ma.masked_all(shape(field))
	length_section	= len(field[0]) - index_break

	#Shift the first part
	new_field[:, :length_section] = field[:, index_break:]

	#Shift the last part
	new_field[:, length_section:] = field[:, :index_break] 

	return new_field

def PeriodicBoundaries3D(lon, lat, field, lon_grids = 1):
    """Add periodic zonal boundaries for 2D field"""
    
    #Empty field with additional zonal boundaries
    lon_2                   = np.zeros((len(lat), len(lon[0]) + lon_grids * 2))
    lat_2                   = np.zeros((len(lat), len(lon_2[0])))
    field_2                 = ma.masked_all((len(field), len(lat), len(lon_2[0])))
    
    #Get the left boundary, which is the right boundary of the original field
    lon_2[:, :lon_grids]       	= lon[:, -lon_grids:]
    lat_2[:, :lon_grids]       	= lat[:, -lon_grids:]
    field_2[:, :, :lon_grids]	= field[:, :, -lon_grids:]
    
    #Same for the right boundary
    lon_2[:, -lon_grids:]        	= lon[:, :lon_grids]
    lat_2[:, -lon_grids:]        	= lat[:, :lon_grids]
    field_2[:, :, -lon_grids:]	= field[:, :, :lon_grids]
    
    #And the complete field
    lon_2[:, lon_grids:-lon_grids]          = lon
    lat_2[:, lon_grids:-lon_grids]          = lat
    field_2[:, :, lon_grids:-lon_grids]     = field
    
    return lon_2, lat_2, field_2

def LowCESMPlot(lon, lat, field):
	"""Returns 4 array's to plot on a global projection"""

	#Left of pole
	lon[lon > 180]	= lon[lon > 180] - 360.0

	lon_1		= lon[:, :160]
	lat_1		= lat[:, :160]
	field_1		= field[:, :160]

	#Right of pole
	lon_2		= lon[:, 159:]
	lat_2		= lat[:, 159:]
	field_2		= field[:, 159:]

	lat_3		= ma.masked_where(lon_2 > 0.0, lat_2)
	field_3		= ma.masked_where(lon_2 > 0.0, field_2)
	lon_3		= ma.masked_where(lon_2 > 0.0, lon_2)

	lon_2[lon_2 < -160] 	= lon_2[lon_2 < -160] + 360
	lat_2			= ma.masked_where(lon_2 < 0.0, lat_2)
	field_2			= ma.masked_where(lon_2 < 0.0, field_2)
	lon_2			= ma.masked_where(lon_2 < 0.0, lon_2)

	#To match at 40W
	index_1		= (fabs(lon[40] - 0.0)).argmin()

	lon_4		= ConverterField(index_1, lon)
	lat_4		= ConverterField(index_1, lat)
	field_4		= ConverterField(index_1, field)

	lon_4		= lon_4[:, 280:300]
	lat_4		= lat_4[:, 280:300]
	field_4		= field_4[:, 280:300]

	return lon_1, lat_1, field_1, lon_2, lat_2, field_2, lon_3, lat_3, field_3, lon_4, lat_4, field_4

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

fh 	= netcdf.Dataset(directory+'Data/Arctic_sea_ice_fraction_March/CESM_year_0001-0050.nc', 'r')

lon		= fh.variables['lon'][:]
lat		= fh.variables['lat'][:] 			
fraction_1	= np.mean(fh.variables['Fraction'][:], axis = 0)

fh.close()

fh 	= netcdf.Dataset(directory+'Data/Arctic_sea_ice_fraction_March/CESM_year_2151-2200.nc', 'r')

fraction_2	= np.mean(fh.variables['Fraction'][:], axis = 0)

fh.close()

#-----------------------------------------------------------------------------------------

lon_1, lat_1, fraction_1_plot, lon_2, lat_2, fraction_2_plot, lon_3, lat_3, fraction_3_plot, lon_4, lat_4, fraction_4_plot	= LowCESMPlot(lon, lat, fraction_2)

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Orthographic(0, 90)})

CS      = ax.contourf(lon_1, lat_1, fraction_1_plot, levels = np.arange(15, 100.01, 4.25), cmap = 'tab20b', transform=ccrs.PlateCarree())
CS      = ax.contourf(lon_2, lat_2, fraction_2_plot, levels = np.arange(15, 100.01, 4.25), cmap = 'tab20b', transform=ccrs.PlateCarree())
CS      = ax.contourf(lon_3, lat_3, fraction_3_plot, levels = np.arange(15, 100.01, 4.25), cmap = 'tab20b', transform=ccrs.PlateCarree())
CS      = ax.contourf(lon_4, lat_4, fraction_4_plot, levels = np.arange(15, 100.01, 4.25), cmap = 'tab20b', transform=ccrs.PlateCarree())

lon_1, lat_1, fraction_1_plot, lon_2, lat_2, fraction_2_plot, lon_3, lat_3, fraction_3_plot, lon_4, lat_4, fraction_4_plot	= LowCESMPlot(lon, lat, fraction_1)

CS_1	= ax.contour(lon_1, lat_1, fraction_1_plot, levels = [15], colors = 'darkblue', linewidths = 2, transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon_2, lat_2, fraction_2_plot, levels = [15], colors = 'darkblue', linewidths = 2, transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon_3, lat_3, fraction_3_plot, levels = [15], colors = 'darkblue', linewidths = 2, transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon_4, lat_4, fraction_4_plot, levels = [15], colors = 'darkblue', linewidths = 2, transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = [15, 40, 60, 80, 100], cax=ax_cb)
cbar.set_label('Sea-ice fraction ($\%$)', fontsize = 12)

ax.gridlines(zorder=10)
ax.add_feature(cfeature.LAND, zorder=10)
ax.coastlines()

ax.scatter([-90, 0, 90, 180], [35, 35, 35, 35], marker = 'o', edgecolor = 'w' , s = 6, facecolors='none', transform=ccrs.PlateCarree())

graphs	= ax.plot([50, 50], [90, 90], '-', color = 'darkblue', linewidth = 2.0, label = '1 - 50', transform=ccrs.PlateCarree(), zorder =0)

legend_labels = [l.get_label() for l in graphs]

ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1).set_zorder(12)


ax.set_title('a) Arctic sea-ice extent, March (2151 - 2200)')

show()



