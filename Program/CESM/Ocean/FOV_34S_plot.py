#Program plots the FOV

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
from scipy.interpolate import CubicSpline

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

def ReadinData(filename):

    fh = netcdf.Dataset(filename, 'r')
    
    time            = fh.variables['time'][:]
    FOV             = fh.variables['F_OV'][:]       #Fresh water
    
    fh.close()
    
    return time, FOV

def ReadinDataGyre(filename):

    fh = netcdf.Dataset(filename, 'r')
    
    time            = fh.variables['time'][:]
    FOV             = fh.variables['F_gyre'][:]     #Fresh water
    
    fh.close()
    
    return time, FOV

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

time, FOV_34S   = ReadinData(directory+'Ocean/FOV_index_section_34S.nc')

#----------------------------------------------------------------------------------------

year_average       = 50
time_average_CS    = int(len(time) / year_average)
time_average_CS    = ma.masked_all((year_average, time_average_CS))
FOV_average_CS     = ma.masked_all((year_average, len(time_average_CS[0])))


for year_start in range(year_average):
    for time_i in range(len(time_average_CS[0])):
        #Loop over each time slice
        if len(time[year_start+time_i*year_average:year_start+(time_i+1)*year_average]) != year_average:
            #No complete period
            continue
            
        #Determine the averages for the cubic splines knots 
        time_average_CS[year_start, time_i] = np.mean(time[year_start+time_i*year_average:year_start+(time_i+1)*year_average])    
        FOV_average_CS[year_start, time_i]  = np.mean(FOV_34S[year_start+time_i*year_average:year_start+(time_i+1)*year_average])    


FOV_mean     = np.zeros((year_average, len(time)))

for year_start in range(year_average):
	#Get the different knots
	time_knot       = time_average_CS[year_start]
	FOV_knot        = FOV_average_CS[year_start]

	#Get non-masked data
	mask_index      = np.where(time_knot.mask == False)[0]
	time_knot       = time_knot[mask_index]
	FOV_knot        = FOV_knot[mask_index]

	#Fit the cubic splines
	cs                   = CubicSpline(time_knot, FOV_knot, bc_type = 'natural')
	FOV_mean[year_start] = cs(time)

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-100, 2500], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')
ax.plot([0.1/0.0003, 0.1/0.0003], [-0.4, 0.312], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.2/0.0003, 0.2/0.0003], [-0.4, 0.312], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.3/0.0003, 0.3/0.0003], [-0.4, 0.312], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.4/0.0003, 0.4/0.0003], [-0.4, 0.312], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.5/0.0003, 0.5/0.0003], [-0.4, 0.312], linestyle = '--', color = 'c', linewidth = 1)
ax.plot([0.6/0.0003, 0.6/0.0003], [-0.4, 0.312], linestyle = '--', color = 'c', linewidth = 1)

ax.text(0.1/0.0003, 0.312, '0.1 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.2/0.0003, 0.312, '0.2 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.3/0.0003, 0.312, '0.3 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.4/0.0003, 0.312, '0.4 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.5/0.0003, 0.312, '0.5 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)
ax.text(0.6/0.0003, 0.312, '0.6 Sv', verticalalignment='bottom', horizontalalignment='center', color = 'c', fontsize=11)

ax.plot(time, FOV_34S, '-k', linewidth = 0.5)

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(1, 2200)
ax.set_ylim(-0.35, 0.35)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.20, 0.17, 0.30, 0.28])


graph_34S	= ax2.plot(time, FOV_34S, '-k', linewidth = 1, label = '$F_{\mathrm{ovS}}$')
graph_mean	= ax2.plot(time, np.mean(FOV_mean, axis = 0), '-r', linewidth = 2, label = '$F_{\mathrm{ovS}}$ mean')

#-----------------------------------------------------------------------------------------

ax2.set_xlim(1650, 1850)
ax2.set_ylim(-0.2, 0)
ax2.grid()

ax2.set_yticks([-0.2, -0.1, 0])

ax2.set_title('$F_{\mathrm{ovS}}$ minimum', fontsize = 10.5)

graphs	      	= graph_34S + graph_mean
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax2.legend(graphs, legend_labels, loc=(0.8, 0.02), ncol=1, framealpha = 1.0, numpoints = 1, prop={'size': 9})

ax2.quiver(1732.48, -0.05, 0, -1.2, scale = 5, color = 'r', zorder = 10, width = 0.01)
ax2.text(1732.48, -0.04, '$\mathrm{d}_t F_{\mathrm{ovS}} = 0$', verticalalignment='bottom', horizontalalignment='center', color = 'r', fontsize = 9)

#-----------------------------------------------------------------------------------------

ax.set_title('a) $F_{\mathrm{ovS}}$')

show()


