#Program plots the freshwater transport and AMOC strength

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from matplotlib.collections import LineCollection

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	FOV		    = fh.variables['F_OV'][:]	    #Freshwater (Sv)

	fh.close()

	return time, FOV

def ReadinDataAMOC(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	AMOC		= fh.variables['Transport'][:]	#AMOC strength

	fh.close()

	return time, AMOC

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

depth_min 	= 0
depth_max	= 1000


time, FOV	= ReadinData(directory+'Ocean/FOV_index_section_34S.nc')
time, AMOC	= ReadinDataAMOC(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')

year_average 	= 25

time_average    	= int(len(time) / year_average)
FOV_average     	= np.zeros(time_average)
AMOC_average    	= np.zeros(time_average)

for time_i in range(time_average):
	#Loop over each time slice
	FOV_average[time_i] 		= np.mean(FOV[time_i*year_average:(time_i+1)*year_average])
	AMOC_average[time_i] 		= np.mean(AMOC[time_i*year_average:(time_i+1)*year_average])

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-0.28, -0.05], y1 = -100, y2 = 100, alpha=0.25, edgecolor='orange', facecolor='orange')
ax.fill_between([-1, 1], y1 = 16, y2 = 19, alpha=0.25, edgecolor='orange', facecolor='orange')

graph_2 = ax.plot([-0.24464943, 0.49348203], [11.26461723, 22.5834036], '-c', linewidth = 1.5, zorder = 0, label = 'CMIP6')
ax.plot([-0.24464943, 0.49348203], [11.26461723+2.704060858924384, 22.5834036+2.704060858924384], '--c')
ax.plot([-0.24464943, 0.49348203], [11.26461723-2.704060858924384, 22.5834036-2.704060858924384], '--c')

graph_1 = ax.plot(FOV_average[:88], AMOC_average[:88], '-k', linewidth = 1.5, label = 'CESM')
ax.scatter(0.5*(FOV_average[0]+FOV_average[1]), 0.5*(AMOC_average[0]+AMOC_average[1]), s = 50, color = 'k', marker = 'o', edgecolor = 'gray', zorder = 10)
ax.scatter(0.5*(FOV_average[40]+FOV_average[41]), 0.5*(AMOC_average[40]+AMOC_average[41])-0.1, s = 60, color = 'k', marker = 'v', edgecolor = 'gray', zorder = 10)
ax.scatter(0.5*(FOV_average[68]+FOV_average[69]), 0.5*(AMOC_average[68]+AMOC_average[69]), s = 50, color = 'k', marker = 'D', edgecolor = 'gray', zorder = 10)
ax.scatter(0.5*(FOV_average[86]+FOV_average[87]), 0.5*(AMOC_average[86]+AMOC_average[87]), s = 50, color = 'k', marker = 's', edgecolor = 'gray', zorder = 10)

ax.text(0.5*(FOV_average[0]+FOV_average[1])+0.00, 0.5*(AMOC_average[0]+AMOC_average[1])-1.2, '1-50', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 10)
ax.text(0.5*(FOV_average[40]+FOV_average[41]), 0.5*(AMOC_average[40]+AMOC_average[41])-1.3, '1001-1050', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 10)
ax.text(0.5*(FOV_average[68]+FOV_average[69])+0.065, 0.5*(AMOC_average[68]+AMOC_average[69]), '1701-1750', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 10)
ax.text(0.5*(FOV_average[86]+FOV_average[87])+0.00, 0.5*(AMOC_average[86]+AMOC_average[87])-1.2, '2151-2200', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 10)

ax.set_xlim(-0.35, 0.35)
ax.set_ylim(-2, 22)
ax.set_xlabel('$F_{\mathrm{ovS}}$ (Sv)')
ax.set_ylabel('AMOC strength (Sv)')
ax.grid()

graphs	      = graph_1 + graph_2

legend_labels = [l.get_label() for l in graphs]

ax.legend(graphs, legend_labels, loc='upper left', fancybox=True, scatterpoints=1, ncol = 1, framealpha = 1.0)
ax.set_title('d) $F_{\mathrm{ovS}}$ and AMOC strength')

show()