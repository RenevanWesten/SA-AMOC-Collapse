#Program plots the freshwater transport and the various components

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

def ReadinDataFOV(filename):

	fh = netcdf.Dataset(filename, 'r')

	time	= fh.variables['time'][:]		
	FOV		= fh.variables['F_OV'][:]	#Freshwater overturning component

	fh.close()

	return time, FOV

def ReadinDataGyre(filename):

	fh = netcdf.Dataset(filename, 'r')

	time	= fh.variables['time'][:]		
	FOV		= fh.variables['F_gyre'][:]	#Freshwater transport by gyre

	fh.close()

	return time, FOV

def ReadinDataFW(filename):

	fh = netcdf.Dataset(filename, 'r')

	time	= fh.variables['time'][:2200]		
	FOV		= fh.variables['FW_transport'][:2200]	#Freshwater transport

	fh.close()

	return time, FOV
    
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

time, FW_34S		= ReadinDataFW(directory+'Ocean/Freshwater_transport_section_34S.nc')
time, FOV_34S		= ReadinDataFOV(directory+'Ocean/FOV_index_section_34S.nc')
time, FW_gyre_34S 	= ReadinDataGyre(directory+'Ocean/FW_gyre_section_34S.nc')

time, FW_60N		= ReadinDataFW(directory+'Ocean/Freshwater_transport_section_60N.nc')
time, FOV_60N		= ReadinDataFOV(directory+'Ocean/FOV_index_section_60N.nc')
time, FW_gyre_60N 	= ReadinDataGyre(directory+'Ocean/FW_gyre_section_60N.nc')

time, FW_med 		= ReadinDataFW(directory+'Ocean/Freshwater_transport_Mediterranean.nc')

#-----------------------------------------------------------------------------------------

#Determine the freshwater transport into the Atlantic Ocean
#A positive transport means net freshwater through the section, thus 60N and the Mediterranean are subtracted
FW_transport	= FW_34S - FW_60N - FW_med

FW_transport_diff	= np.mean(FW_transport[1700:1750])-np.mean(FW_transport[0:50])
FOV_34S_diff		= np.mean(FOV_34S[1700:1750])-np.mean(FOV_34S[0:50])
FOV_60N_diff		= -(np.mean(FOV_60N[1700:1750])-np.mean(FOV_60N[0:50]))
FW_gyre_34S_diff	= (np.mean(FW_gyre_34S[1700:1750])-np.mean(FW_gyre_34S[0:50]))
FW_gyre_60N_diff	= -(np.mean(FW_gyre_60N[1700:1750])-np.mean(FW_gyre_60N[0:50]))

print(FW_transport_diff)
print('FOV (34S) contribution:', FOV_34S_diff / FW_transport_diff*100.)
print('FOV (60N) contribution:', FOV_60N_diff / FW_transport_diff*100.)
print('Gyre (34S) contribution:', FW_gyre_34S_diff / FW_transport_diff * 100.)
print('Gyre (60N) contribution:', FW_gyre_60N_diff / FW_transport_diff * 100.)

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(time, FOV_60N, '-', color = 'gray', linewidth = 0.5, label = '$F_{\mathrm{ovN}}$')
ax.plot(time, FW_gyre_60N, '-', color = 'cyan', linewidth = 0.5, label = '$F_{\mathrm{azN}}$')
ax.plot(time, FW_60N, '-', color = 'firebrick', linewidth = 0.5, label = r'$F_{\nabla\mathrm{N}}$')

ax.plot(time, FW_34S, '-r', linewidth = 0.5, label = r'$F_{\nabla\mathrm{S}}$')
ax.plot(time, FW_gyre_34S, '-b', linewidth = 0.5, label = '$F_{\mathrm{azS}}$')
ax.plot(time, FOV_34S, '-k', linewidth = 0.5, label = '$F_{\mathrm{ovS}}$')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(1, 2200)
ax.set_ylim(-0.55, 0.55)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])

graph_1	= ax.plot(time, FOV_34S-100, '-k', linewidth = 2, label = '$F_{\mathrm{ovS}}$')
graph_2	= ax.plot(time, FW_gyre_34S-100, '-b', linewidth = 2, label = '$F_{\mathrm{azS}}$')
graph_3	= ax.plot(time, FW_34S-100, '-r', linewidth = 2, label = r'$F_{\nabla\mathrm{S}}$')

graph_4	= ax.plot(time, FOV_60N-100, '-', color = 'gray', linewidth = 2, label = '$F_{\mathrm{ovN}}$')
graph_5	= ax.plot(time, FW_gyre_60N-100, '-', color = 'cyan', linewidth = 2, label = '$F_{\mathrm{azN}}$')
graph_6	= ax.plot(time, FW_60N-100, '-', color = 'firebrick', linewidth = 2, label = r'$F_{\nabla\mathrm{N}}$')


graphs	      = graph_1 + graph_2 + graph_3

legend_labels = [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

graphs	      = graph_4 + graph_5 + graph_6

legend_labels = [l.get_label() for l in graphs]
legend_2      = ax.legend(graphs, legend_labels, loc='lower right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.add_artist(legend_1)

#-----------------------------------------------------------------------------------------

ax.set_title('c) Freshwater transport at boundaries')

show()



