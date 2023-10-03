#Program computes the meridional heat transport

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy.interpolate import interp1d

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

def ReadinDataAtmosphere(filename):

	fh = netcdf.Dataset(filename, 'r')

	lat	= fh.variables['lat'][:]		#Latitudes (degrees N)
	MHT	= fh.variables['MHT'][:]		#Meridional heat transport (PW)

	fh.close()

	return lat, MHT


def ReadinDataOcean(filename):

	fh = netcdf.Dataset(filename, 'r')

	lat = fh.variables['lat'][:]		                        #Latitudes (degrees N)
	MHT = np.mean(fh.variables['MHT_Global'][:], axis = 0)		#Meridional heat flux (PW)

	fh.close()

	return lat, MHT
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

lat_atm, MHT_all_1		= ReadinDataAtmosphere(directory+'Atmosphere/Meridional_heat_transport_year_0001-0050.nc')
lat_atm, MHT_all_2		= ReadinDataAtmosphere(directory+'Atmosphere/Meridional_heat_transport_year_2151-2200.nc')
lat_ocn, MHT_ocn_1		= ReadinDataOcean(directory+'Data/AMOC_MHT/CESM_year_0001-0050.nc')
lat_ocn, MHT_ocn_2		= ReadinDataOcean(directory+'Data/AMOC_MHT/CESM_year_2151-2200.nc')
#-----------------------------------------------------------------------------------------

MHT_ocn_1		= interp1d(lat_ocn, MHT_ocn_1)(np.array(lat_atm[8:-8]))
MHT_ocn_2		= interp1d(lat_ocn, MHT_ocn_2)(np.array(lat_atm[8:-8]))

MHT_atm_1	    = np.copy(MHT_all_1)
MHT_atm_2	    = np.copy(MHT_all_2)
MHT_atm_1[8:-8] -= MHT_ocn_1
MHT_atm_2[8:-8] -= MHT_ocn_2
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_1		= ax.plot(lat_atm, MHT_all_1, '-r', linewidth = 1.5, label = '1-50')
graph_2		= ax.plot(lat_atm, MHT_all_2, '-b', linewidth = 1.5, label = '2151-2200')

ax.plot(lat_atm, MHT_atm_1, '--r', linewidth = 1.5, label = '1-50')
ax.plot(lat_atm, MHT_atm_2, '--b', linewidth = 1.5, label = '2151-2200')

ax.set_ylabel('Meridional heat transport (PW)')
ax.set_xlim(-70, 70)
ax.set_ylim(-6, 6)
ax.grid()

ax.set_xticks(np.arange(-60, 60.1, 20))
ax.set_xticklabels(['60$^{\circ}$S', '40$^{\circ}$S', '20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N'])
ax.set_yticks(np.arange(-6, 6.1, 2))

ax2 	= ax.twinx()

graph_3	= ax2.plot(lat_atm, MHT_all_2 - MHT_all_1, '-k', linewidth = 2.0, label = '$\Delta$MHT')
ax2.plot(lat_atm, MHT_atm_2 - MHT_atm_1, '--k', linewidth = 1.5, label = '$\Delta$MHT')

ax2.set_ylim(-1, 1)
ax2.set_ylabel('Meridional heat transport difference (PW)')
ax2.set_yticks([-1, -0.5, 0, 0.5, 1])

graphs	      	= graph_1 + graph_2 + graph_3

legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='lower right', ncol=1, framealpha=1, numpoints = 1)

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = 'Total')
graph_2		= ax.plot([-100, -100], [-100, -100], '--k', linewidth = 1.5, label = 'Atmosphere')

graphs	      	= graph_1 + graph_2

legend_labels 	= [l.get_label() for l in graphs]
legend_2	= ax.legend(graphs, legend_labels, loc = 'upper left', ncol=1, framealpha = 1.0)
ax.add_artist(legend_1)

ax.set_title('b) Total and atmospheric meridional heat transport')

show()

