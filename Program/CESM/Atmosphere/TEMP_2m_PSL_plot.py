#Program plots the 2-meter temperature and sea-level pressure difference

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

fh 	= netcdf.Dataset(directory+'Atmosphere/TEMP_2m_PSL_year_0001-0050.nc', 'r')

lon	    = fh.variables['lon'][:]
lat	    = fh.variables['lat'][:] 			
temp_1	= fh.variables['TEMP_2m'][:] 
pres_1	= np.mean(fh.variables['PSL'][:], axis = 0)

fh.close()

fh 	= netcdf.Dataset(directory+'Atmosphere/TEMP_2m_PSL_year_2151-2200.nc', 'r')
			
temp_2	= fh.variables['TEMP_2m'][:] 
pres_2	= np.mean(fh.variables['PSL'][:], axis = 0)

fh.close()

#-----------------------------------------------------------------------------------------
#Rescale the temperature plot
scale	= 5.0
cut_off	= 5

temp_plot			            = np.mean(temp_2, axis = 0) - np.mean(temp_1, axis = 0)
temp_plot[temp_plot < -cut_off]	= (temp_plot[temp_plot < -cut_off] - -cut_off) / scale - cut_off
temp_plot[temp_plot > cut_off]	= (temp_plot[temp_plot > cut_off] - cut_off) / scale + cut_off
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Robinson()})

CS      = ax.contourf(lon, lat, temp_plot, levels = np.arange(-8, 8.01, 0.5), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = [-8, -6, -4, -2, 0, 2, 4, 6, 8], cax=ax_cb)
cbar.ax.set_yticklabels([-20, -10, -4, -2, 0, 2, 4, 10, 20])
cbar.set_label('2-meter temperature difference ($^{\circ}$C)')

CS_1	= ax.contour(lon, lat, pres_2 - pres_1, levels = [0], colors = 'k', linewidths = 2, transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon, lat, pres_2 - pres_1, levels = [-1], colors = 'c', linewidths = 2, transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon, lat, -(pres_2 - pres_1), levels = [2], colors = 'b', linewidths = 2, transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon, lat, -(pres_2 - pres_1), levels = [-1], colors = 'firebrick', linewidths = 2, transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon, lat, pres_2 - pres_1, levels = [2], colors = 'r', linewidths = 2, transform=ccrs.PlateCarree())

ax.set_global()
ax.gridlines()
ax.coastlines()

for lat_i in range(0, len(lat), 3):
	for lon_i in range(0, len(lon), 3):
		#Determine significant difference
		p_value = Welch(temp_1[:, lat_i, lon_i], temp_2[:, lat_i, lon_i])

		if p_value <= 0.95:
			#Non-significant difference
			ax.scatter(lon[lon_i], lat[lat_i], marker = 'o', edgecolor = 'k' , s = 6, facecolors='none', transform=ccrs.PlateCarree())

ax.set_title('a) 2-meter temperature and sea-level pressure')
show()

