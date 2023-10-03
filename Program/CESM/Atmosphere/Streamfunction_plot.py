#Program plots the atmospheric streamfunction difference

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
#-----------------------------------------------------------------------------------------
fh 	= netcdf.Dataset(directory+'Atmosphere/Streamfunction_year_0001-0050.nc', 'r')

pres		= fh.variables['lev'][:]
lat		    = fh.variables['lat'][:]				
stream_1	= fh.variables['SF'][:] / 10**10.0	#Streamfunction (10^10 kg / s)	

fh.close()

fh 	= netcdf.Dataset(directory+'Atmosphere/Streamfunction_year_2151-2200.nc', 'r')
			
stream_2	= fh.variables['SF'][:] / 10**10.0	#Streamfunction (10^10 kg / s)	

fh.close()

#-----------------------------------------------------------------------------------------

#Take the ensemble mean
stream_plot		= np.mean(stream_2, axis = 0) - np.mean(stream_1, axis = 0)
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

CS	= ax.contourf(lat, pres, stream_plot, levels = np.arange(-4, 4.1, 0.25), extend = 'both', cmap = 'PuOr_r')
cbar	= colorbar(CS, ticks = np.arange(-4, 4.1, 1))
cbar.set_label('Mass streamfunction difference ($10^{10}$ kg s$^{-1}$)')

ax.set_ylabel('Pressure (hPa)')
ax.set_xlim(-70, 70)
ax.set_ylim(950, 50)
ax.set_yscale('log')

ax.set_xticks(np.arange(-60, 60.1, 20))
ax.set_xticklabels(['60$^{\circ}$S', '40$^{\circ}$S', '20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])
ax.set_yticks([50, 100, 200, 300, 400, 500, 600, 850])
ax.set_yticklabels([50, 100, 200, 300, '', 500, '', 850])

CS3	= ax.contour(lat, pres, np.mean(stream_1, axis = 0) - 2, levels = [-1], colors = 'firebrick', linewidths = 2)
CS3	= ax.contour(lat, pres, np.mean(stream_1, axis = 0), levels = [2], colors = 'r', linewidths = 2)
CS3	= ax.contour(lat, pres, np.mean(stream_1, axis = 0), levels = [-1], colors = 'cyan', linewidths = 2)
CS3	= ax.contour(lat, pres, np.mean(stream_1, axis = 0) + 4, levels = [2], colors = 'b', linewidths = 2)

for lev_i in range(len(pres)):
	for lat_i in range(len(lat)):
		#Determine significant difference

		p_value = Welch(stream_1[:, lev_i, lat_i], stream_2[:, lev_i, lat_i])

		if p_value <= 0.95:
			#Non-significant difference
			ax.scatter(lat[lat_i], pres[lev_i], marker = 'o', edgecolor = 'k' , s = 6, facecolors='none', zorder = 10)


ax.set_title('e) Meridional mass streamfunction')
show()
