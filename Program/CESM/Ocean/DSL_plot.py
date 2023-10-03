#Program plots the DSL difference

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

	lat_2		= ma.masked_where(lon_2 < 0.0, lat_2)
	field_2		= ma.masked_where(lon_2 < 0.0, field_2)
	lon_2		= ma.masked_where(lon_2 < 0.0, lon_2)

	#To match at 40W
	index_1		= (fabs(lon[40] - 0.0)).argmin()

	lon_4		= ConverterField(index_1, lon)
	lat_4		= ConverterField(index_1, lat)
	field_4		= ConverterField(index_1, field)

	lon_4		= lon_4[:, 280:300]
	lat_4		= lat_4[:, 280:300]
	field_4		= field_4[:, 280:300]

	return lon_1, lat_1, field_1, lon_2, lat_2, field_2, lon_3, lat_3, field_3, lon_4, lat_4, field_4

def Welch(data_1, data_2):
	"""Conducts the Welch t-test"""
	
	#Determine the means
	mean_1	= np.mean(data_1)
	mean_2	= np.mean(data_2)
	
	#Determine the corrected sample standard deviations
	std_1	= np.sqrt(1.0 / (len(data_1) - 1) * np.sum((data_1 - mean_1)**2.0))
	std_2	= np.sqrt(1.0 / (len(data_2) - 1) * np.sum((data_2 - mean_2)**2.0))

	#Determine the Welch t-value
	t_welch	= (mean_1 - mean_2) / np.sqrt((std_1**2.0 / len(data_1)) + (std_2**2.0 / len(data_2)))

	#Determine the degrees of freedome (dof)
	dof	= ((std_1**2.0 / len(data_1)) + (std_2**2.0 / len(data_2)))**2.0 / ((std_1**4.0 / (len(data_1)**2.0 * (len(data_1) - 1))) + (std_2**4.0 / (len(data_2)**2.0 * (len(data_2) - 1))))

	#Get the significance levels and the corresponding critical values (two-sided)
	sig_levels 	= np.arange(50, 100, 0.5) / 100.0
	t_crit		= stats.t.ppf((1.0 + sig_levels) / 2.0, dof)

	#Get the indices where the significance is exceeding the critical values
	sig_index	= np.where(fabs(t_welch) > t_crit)[0]
	significant	= 0.0

	if len(sig_index) > 0:
		#If there are significance values, take the highest significant level
		significant = sig_levels[sig_index[-1]]

	return significant


#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

fh 	= netcdf.Dataset(directory+'Ocean/DSL_fields.nc', 'r')

lon	    = fh.variables['lon'][:]
lat	    = fh.variables['lat'][:] 			
DSL_1	= fh.variables['DSL_0001-0050'][:] 
DSL_2	= fh.variables['DSL_2151-2200'][:] 
p_value = fh.variables['DSL_diff'][:]

fh.close()

#-----------------------------------------------------------------------------------------
#Rescale the dyanmic sea-level plot
scale	= 2
cut_off	= 20

DSL_plot			            = DSL_2 - DSL_1
DSL_plot[DSL_plot < -cut_off]	= (DSL_plot[DSL_plot < -cut_off] - -cut_off) / scale - cut_off
DSL_plot[DSL_plot > cut_off]	= (DSL_plot[DSL_plot > cut_off] - cut_off) / scale + cut_off
#-----------------------------------------------------------------------------------------

lon_1, lat_1, DSL_1_plot, lon_2, lat_2, DSL_2_plot, lon_3, lat_3, DSL_3_plot, lon_4, lat_4, DSL_4_plot	= LowCESMPlot(lon, lat, DSL_plot)

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Robinson()})

CS      = ax.contourf(lon_1, lat_1, DSL_1_plot, levels = np.arange(-60, 60.01, 2.5), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())
CS      = ax.contourf(lon_2, lat_2, DSL_2_plot, levels = np.arange(-60, 60.01, 2.5), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())
CS      = ax.contourf(lon_3, lat_3, DSL_3_plot, levels = np.arange(-60, 60.01, 2.5), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())
CS      = ax.contourf(lon_4, lat_4, DSL_4_plot, levels = np.arange(-60, 60.01, 2.5), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())


divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = [-60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60], cax=ax_cb)
cbar.ax.set_yticklabels([-100, -80, -60, -40, -20, -10, 0, 10, 20, 40, 60, 80, 100])
cbar.set_label('Dynamic sea-level difference (cm)')

ax.set_global()
ax.gridlines(zorder = 11)
ax.add_feature(cfeature.LAND, zorder=10)
ax.coastlines()

for lat_i in range(0, len(lat), 10):
	for lon_i in range(0, len(lon[0]), 10):
		#Determine significant difference

		if DSL_plot[lat_i, lon_i] is ma.masked:
			continue

		if p_value[lat_i, lon_i] <= 0.95:
			#Non-significant difference
			ax.scatter(lon[lat_i, lon_i], lat[lat_i, lon_i], marker = 'o', edgecolor = 'k' , s = 6, facecolors='none', transform=ccrs.PlateCarree())

ax.set_title('d) Dynamic sea level')
show()

