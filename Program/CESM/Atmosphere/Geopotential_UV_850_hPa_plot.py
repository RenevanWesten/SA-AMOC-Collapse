#Program plots the 850 hPa geopoential and horizontal velocity difference

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

fh 	= netcdf.Dataset(directory+'Atmosphere/Z3_UV_850_hPa_year_0001-0050.nc', 'r')

lon	= fh.variables['lon'][:]
lat	= fh.variables['lat'][:] 			
Z3_1	= fh.variables['Z3_850'][:] 
u_vel_1	= fh.variables['U_850'][:]
v_vel_1	= fh.variables['V_850'][:]

fh.close()

fh 	= netcdf.Dataset(directory+'Atmosphere/Z3_UV_850_hPa_year_2151-2200.nc', 'r')
			
Z3_2	= fh.variables['Z3_850'][:] 
u_vel_2	= fh.variables['U_850'][:]
v_vel_2	= fh.variables['V_850'][:]

fh.close()

#-----------------------------------------------------------------------------------------
#Rescale the precerature plot
scale	= 5.0
cut_off	= 5

Z3_plot		= np.mean(Z3_2, axis = 0) - np.mean(Z3_1, axis = 0)
u_vel_plot	= np.mean(u_vel_2, axis = 0) - np.mean(u_vel_1, axis = 0)
v_vel_plot	= np.mean(v_vel_2, axis = 0) - np.mean(v_vel_1, axis = 0)

#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Robinson()})

CS      = ax.contourf(lon, lat, Z3_plot, levels = np.arange(-30, 30.01, 2), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-30, 30.01, 10), cax=ax_cb)
cbar.set_label('Geopotential height difference (m)')

ax.set_global()
ax.add_feature(cfeature.LAND)
ax.gridlines()
ax.coastlines()

for lat_i in range(0, len(lat), 3):
	for lon_i in range(0, len(lon), 3):
		#Determine significant difference

		if Z3_plot[lat_i, lon_i] is ma.masked:
			continue

		p_value = Welch(Z3_1[:, lat_i, lon_i], Z3_2[:, lat_i, lon_i])

		if p_value <= 0.95:
			#Non-significant difference
			ax.scatter(lon[lon_i], lat[lat_i], marker = 'o', edgecolor = 'k' , s = 6, facecolors='none', transform=ccrs.PlateCarree())


scale_arrow	= 4
Q = ax.quiver(lon[::scale_arrow], lat[::scale_arrow], u_vel_plot[::scale_arrow, ::scale_arrow], v_vel_plot[::scale_arrow, ::scale_arrow], scale = 40, transform=ccrs.PlateCarree())

qk = ax.quiverkey(Q, 0.17, 0.28, 2, '2 m s$^{-1}$', labelpos = 'S', coordinates='figure')

ax.set_title('c) 850 hPa geopotential height and horizontal velocity')
show()

