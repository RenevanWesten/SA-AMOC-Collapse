    #Program plots the oceanic meridional heat transport

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	lat		        = fh.variables['lat'][:]		                        #Latitudes (degrees N)
	MHT_atlantic	= np.mean(fh.variables['MHT_Atlantic'][:], axis = 0)	#Meridional heat flux (PW)
	MHT_pacific	    = np.mean(fh.variables['MHT_Indian-Pacific'][:], axis = 0)	#Meridional heat flux (PW)
	MHT_global	    = np.mean(fh.variables['MHT_Global'][:], axis = 0)		#Meridional heat flux (PW)

	fh.close()

	return lat, MHT_atlantic, MHT_pacific, MHT_global

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

lat, MHT_atlantic_1, MHT_pacific_1, MHT_global_1	= ReadinData(directory+'Data/AMOC_MHT/CESM_year_0001-0050.nc')
lat, MHT_atlantic_2, MHT_pacific_2, MHT_global_2	= ReadinData(directory+'Data/AMOC_MHT/CESM_year_2151-2200.nc')

print('Meridional heat transport difference (%) at 26N')
print((MHT_atlantic_2[202] - MHT_atlantic_1[202]) / MHT_atlantic_1[202] * 100)

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_1		= ax.plot(lat, MHT_global_1, '-r', linewidth = 1.5, label = '1-50')
graph_2		= ax.plot(lat, MHT_global_2, '-b', linewidth = 1.5, label = '2151-2200')

ax.plot(lat, MHT_atlantic_1, '--r', linewidth = 1.5, label = '1-50')
ax.plot(lat, MHT_atlantic_2, '--b', linewidth = 1.5, label = '2151-2200')

ax.plot(lat, MHT_pacific_1, ':r', linewidth = 1.5, label = '1-50')
ax.plot(lat, MHT_pacific_2, ':b', linewidth = 1.5, label = '2151-2200')

ax.set_ylabel('Meridional heat transport (PW)')
ax.set_xlim(-70, 70)
ax.set_ylim(-2, 2)
ax.grid()

ax.set_xticks(np.arange(-60, 60.1, 20))
ax.set_xticklabels(['60$^{\circ}$S', '40$^{\circ}$S', '20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])
ax.set_yticks(np.arange(-2, 2.1, 1))

ax2 	= ax.twinx()

graph_3	= ax2.plot(lat, MHT_global_2 - MHT_global_1, '-k', linewidth = 1.5, label = '$\Delta$MHT')
ax2.plot(lat, MHT_atlantic_2 - MHT_atlantic_1, '--k', linewidth = 1.5, label = '$\Delta$MHT')
ax2.plot(lat, MHT_pacific_2 - MHT_pacific_1, ':k', linewidth = 1.5, label = '$\Delta$MHT')

ax2.set_ylim(-1, 1)
ax2.set_ylabel('Meridional heat transport difference (PW)')
ax2.set_yticks([-1, -0.5, 0, 0.5, 1])

graphs	      = graph_1 + graph_2 + graph_3	

legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha=1, numpoints = 1)

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = 'Global')
graph_2		= ax.plot([-100, -100], [-100, -100], '--k', linewidth = 1.5, label = 'Atlantic')
graph_3		= ax.plot([-100, -100], [-100, -100], ':k', linewidth = 1.5, label = 'Indian-Pacific')

graphs	      	= graph_1 + graph_2 + graph_3

legend_labels 	= [l.get_label() for l in graphs]
legend_2	= ax.legend(graphs, legend_labels, loc = 'upper left', ncol=1, framealpha = 1.0)
ax.add_artist(legend_1)

ax.set_title('a) Oceanic meridional heat transport')

show()

