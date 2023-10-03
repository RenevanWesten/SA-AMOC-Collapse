#Program plots the full AMOC and Meridional heat transport

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'/Data/AMOC_structure/CESM_year_0001-0050.nc', 'r')

depth	= fh.variables['depth'][:] 		
lat	    = fh.variables['lat'][:] 		
AMOC_1	= np.mean(fh.variables['AMOC'][:], axis = 0)	

fh.close()

fh = netcdf.Dataset(directory+'/Data/AMOC_structure/CESM_year_1701-1750.nc', 'r')
		
AMOC_2	= np.mean(fh.variables['AMOC'][:], axis = 0)	

fh.close()

fh = netcdf.Dataset(directory+'/Data/AMOC_structure/CESM_year_2151-2200.nc', 'r')

AMOC_3	= np.mean(fh.variables['AMOC'][:], axis = 0)	

fh.close()

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'/Data/AMOC_MHT/CESM_year_0001-0050.nc', 'r')

lat_MHT		    = fh.variables['lat'][:]		                       #Latitudes (degrees N)
MHT_atlantic_1	= np.mean(fh.variables['MHT_Atlantic'][:], axis = 0)	#Meridional heat flux (PW)

fh.close()

fh = netcdf.Dataset(directory+'/Data/AMOC_MHT/CESM_year_1701-1750.nc', 'r')

MHT_atlantic_2	= np.mean(fh.variables['MHT_Atlantic'][:], axis = 0)	#Meridional heat flux (PW)

fh.close()

fh = netcdf.Dataset(directory+'/Data/AMOC_MHT/CESM_year_2151-2200.nc', 'r')

MHT_atlantic_3	= np.mean(fh.variables['MHT_Atlantic'][:], axis = 0)	#Meridional heat flux (PW)

fh.close()

#-----------------------------------------------------------------------------------------


fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1], 'hspace': 0.05})

x, y	= np.meshgrid(lat, depth)
CS	= ax1.contourf(x, y, AMOC_1, levels = np.arange(-6, 21.1, 0.5), extend = 'both', cmap = 'Spectral_r')

fig.subplots_adjust(right=0.8)
cbar 	= fig.add_axes([0.82, 0.278, 0.030, 0.60])
fig.colorbar(CS, cax = cbar, ticks = np.arange(-6, 21.1, 3))
cbar.set_ylabel('Atlantic Meridional Overturning Circulation (Sv)')

CS_1	= ax1.contour(x, y, AMOC_1, levels = [15], colors = 'k', linewidths = 2)
CS_2	= ax1.contour(x, y, AMOC_1, levels = [12], colors = 'r', linewidths = 2)
CS_3	= ax1.contour(x, y, AMOC_1, levels = [9], colors = 'b', linewidths = 2)
CS_3	= ax1.contour(x, y, AMOC_1, levels = [-1], colors = 'k', linewidths = 2)

ax1.set_xlim(-30, 70)
ax1.set_ylim(5200, 0)
ax1.set_ylabel('Depth (m)')
ax1.set_xticklabels([])

CS_1		= ax1.plot([90, 90], [-1, -1], linestyle = '-', color = 'k', linewidth = 2, label = '15 Sv')
CS_2		= ax1.plot([90, 90], [-1, -1], linestyle = '-', color = 'r', linewidth = 2, label = '12 Sv')
CS_3		= ax1.plot([90, 90], [-1, -1], linestyle = '-', color = 'b', linewidth = 2, label = '9 Sv')
CS_4		= ax1.plot([90, 90], [-1, -1], linestyle = '--', color = 'k', linewidth = 2, label = '-1 Sv')

graphs		= CS_1 + CS_2 + CS_3 + CS_4
legend_labels = [l.get_label() for l in graphs]
legend_1      = ax1.legend(graphs, legend_labels, loc = 'lower right', ncol=1, numpoints = 1, framealpha = 1.0)

ax2.plot(lat_MHT, MHT_atlantic_1, '-k', linewidth = 2)

ax2.set_ylabel('MHT (PW)')
ax2.set_xlim(-30, 70)
ax2.set_ylim(-1.2, 1.2)
ax2.grid()

ax2.set_xticks(np.arange(-20, 60.1, 20))
ax2.set_xticklabels(['20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])

ax1.set_title('b) Atlantic overturning circulation and MHT (1 - 50)')
 
#-----------------------------------------------------------------------------------------


fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1], 'hspace': 0.05})

x, y	= np.meshgrid(lat, depth)
CS	= ax1.contourf(x, y, AMOC_2, levels = np.arange(-6, 21.1, 0.5), extend = 'both', cmap = 'Spectral_r')

fig.subplots_adjust(right=0.8)
cbar 	= fig.add_axes([0.82, 0.278, 0.030, 0.60])
fig.colorbar(CS, cax = cbar, ticks = np.arange(-6, 21.1, 3))
cbar.set_ylabel('Atlantic Meridional Overturning Circulation (Sv)')

CS_1	= ax1.contour(x, y, AMOC_2, levels = [15], colors = 'k', linewidths = 2)
CS_2	= ax1.contour(x, y, AMOC_2, levels = [12], colors = 'r', linewidths = 2)
CS_3	= ax1.contour(x, y, AMOC_2, levels = [9], colors = 'b', linewidths = 2)
CS_3	= ax1.contour(x, y, AMOC_2, levels = [-1], colors = 'k', linewidths = 2)

ax1.set_xlim(-30, 70)
ax1.set_ylim(5200, 0)
ax1.set_ylabel('Depth (m)')
ax1.set_xticklabels([])

CS_1		= ax1.plot([90, 90], [-1, -1], linestyle = '-', color = 'k', linewidth = 2, label = '15 Sv')
CS_2		= ax1.plot([90, 90], [-1, -1], linestyle = '-', color = 'r', linewidth = 2, label = '12 Sv')
CS_3		= ax1.plot([90, 90], [-1, -1], linestyle = '-', color = 'b', linewidth = 2, label = '9 Sv')
CS_4		= ax1.plot([90, 90], [-1, -1], linestyle = '--', color = 'k', linewidth = 2, label = '-1 Sv')

graphs		= CS_1 + CS_2 + CS_3 + CS_4
legend_labels = [l.get_label() for l in graphs]
legend_1      = ax1.legend(graphs, legend_labels, loc = 'lower right', ncol=1, numpoints = 1, framealpha = 1.0)

ax2.plot(lat_MHT, MHT_atlantic_2, '-k', linewidth = 2)

ax2.set_ylabel('MHT (PW)')
ax2.set_xlim(-30, 70)
ax2.set_ylim(-1.2, 1.2)
ax2.grid()

ax2.set_xticks(np.arange(-20, 60.1, 20))
ax2.set_xticklabels(['20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])

ax1.set_title('c) Atlantic overturning circulation and MHT (1701 - 1750)')
#-----------------------------------------------------------------------------------------


fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1], 'hspace': 0.05})

x, y	= np.meshgrid(lat, depth)
CS	= ax1.contourf(x, y, AMOC_3, levels = np.arange(-6, 21.1, 0.5), extend = 'both', cmap = 'Spectral_r')

fig.subplots_adjust(right=0.8)
cbar 	= fig.add_axes([0.82, 0.278, 0.030, 0.60])
fig.colorbar(CS, cax = cbar, ticks = np.arange(-6, 21.1, 3))
cbar.set_ylabel('Atlantic Meridional Overturning Circulation (Sv)')

CS_1	= ax1.contour(x, y, AMOC_3, levels = [15], colors = 'k', linewidths = 2)
CS_2	= ax1.contour(x, y, AMOC_3, levels = [12], colors = 'r', linewidths = 2)
CS_3	= ax1.contour(x, y, AMOC_3, levels = [9], colors = 'b', linewidths = 2)
CS_3	= ax1.contour(x, y, AMOC_3, levels = [-1], colors = 'k', linewidths = 2)

ax1.set_xlim(-30, 70)
ax1.set_ylim(5200, 0)
ax1.set_ylabel('Depth (m)')
ax1.set_xticklabels([])

CS_1		= ax1.plot([90, 90], [-1, -1], linestyle = '-', color = 'k', linewidth = 2, label = '15 Sv')
CS_2		= ax1.plot([90, 90], [-1, -1], linestyle = '-', color = 'r', linewidth = 2, label = '12 Sv')
CS_3		= ax1.plot([90, 90], [-1, -1], linestyle = '-', color = 'b', linewidth = 2, label = '9 Sv')
CS_4		= ax1.plot([90, 90], [-1, -1], linestyle = '--', color = 'k', linewidth = 2, label = '-1 Sv')

graphs		= CS_1 + CS_2 + CS_3 + CS_4
legend_labels = [l.get_label() for l in graphs]
legend_1      = ax1.legend(graphs, legend_labels, loc = 'lower right', ncol=1, numpoints = 1, framealpha = 1.0)

ax2.plot(lat_MHT, MHT_atlantic_3, '-k', linewidth = 2)

ax2.set_ylabel('MHT (PW)')
ax2.set_xlim(-30, 70)
ax2.set_ylim(-1.2, 1.2)
ax2.grid()

ax2.set_xticks(np.arange(-20, 60.1, 20))
ax2.set_xticklabels(['20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])

ax1.set_title('d) Atlantic overturning circulation and MHT (2151 - 2200)')
#-----------------------------------------------------------------------------------------

show()