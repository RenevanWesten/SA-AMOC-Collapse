#Program plots the precipitation difference

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

fh 	= netcdf.Dataset(directory+'Atmosphere/PREC_EVAP_year_0001-0050.nc', 'r')

lon	= fh.variables['lon'][:]
lat	= fh.variables['lat'][:] 			
prec_1	= fh.variables['PREC'][:] 

fh.close()

fh 	= netcdf.Dataset(directory+'Atmosphere/PREC_EVAP_year_2151-2200.nc', 'r')
			
prec_2	= fh.variables['PREC'][:] 

fh.close()


#-----------------------------------------------------------------------------------------
prec_plot	= np.mean(prec_2, axis = 0) - np.mean(prec_1, axis = 0)
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Robinson()})

CS      = ax.contourf(lon, lat, prec_plot, levels = np.arange(-2, 2.01, 0.2), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-2, 2.01, 1), cax=ax_cb)
cbar.set_label('Precipitation difference (mm day$^{-1}$)')

ax.set_global()
ax.gridlines()
ax.coastlines()

for lat_i in range(0, len(lat), 3):
	for lon_i in range(0, len(lon), 3):
		#Determine significant difference
		p_value = Welch(prec_1[:, lat_i, lon_i], prec_2[:, lat_i, lon_i])

		if p_value <= 0.95:
			#Non-significant difference
			ax.scatter(lon[lon_i], lat[lat_i], marker = 'o', edgecolor = 'k' , s = 6, facecolors='none', transform=ccrs.PlateCarree())

ax.set_title('d) Precipitation')
show()

