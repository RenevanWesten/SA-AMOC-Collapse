#Program determines the FOV minimum using cubic splines

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats
from scipy.interpolate import CubicSpline

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time	= fh.variables['time'][:]		
	FOV		= fh.variables['F_OV'][:]	#Freshwater transport

	fh.close()

	return time, FOV

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

time, FOV	= ReadinData(directory+'Ocean/FOV_index_section_34S.nc')
              

year_average    = 50
time_average    = int(len(time) / year_average)
time_average    = ma.masked_all((year_average, time_average))
FOV_average     = ma.masked_all((year_average, len(time_average[0])))

for year_start in range(year_average):
    for time_i in range(len(time_average[0])):
        #Loop over each time slice
        if len(time[year_start+time_i*year_average:year_start+(time_i+1)*year_average]) != year_average:
            #No complete period
            continue
            
        time_average[year_start, time_i] = np.mean(time[year_start+time_i*year_average:year_start+(time_i+1)*year_average])    
        FOV_average[year_start, time_i]  = np.mean(FOV[year_start+time_i*year_average:year_start+(time_i+1)*year_average])    

#-----------------------------------------------------------------------------------------	

FOV_mean    = np.zeros((year_average, len(time)))
FOV_der     = np.zeros((year_average, len(time)))
FOV_min     = np.zeros(year_average)

for year_start in range(year_average):
    #Get the different knots
    time_knot       = time_average[year_start]
    FOV_knot        = FOV_average[year_start]

    #Get non-masked data
    mask_index      = np.where(time_knot.mask == False)[0]
    time_knot       = time_knot[mask_index]
    FOV_knot        = FOV_knot[mask_index]
    
    #Fit the cubic splines
    cs                   = CubicSpline(time_knot, FOV_knot, bc_type = 'natural')
    FOV_mean[year_start] = cs(time)
    FOV_der[year_start]  = cs(time, 1)

    #Get the minimum
    FOV_min[year_start]	= time[np.argmin(FOV_mean[year_start])]


print('FOV minimum (mean):', np.mean(FOV_min))
print('FOV minimum (5 percentile):', np.percentile(FOV_min, 5))
print('FOV minimum (10 percentile):', np.percentile(FOV_min, 10))
print('FOV minimum (25 percentile):', np.percentile(FOV_min, 25))
print('FOV minimum (75 percentile):', np.percentile(FOV_min, 75))
print('FOV minimum (90 percentile):', np.percentile(FOV_min, 90))
print('FOV minimum (95 percentile):', np.percentile(FOV_min, 95))

#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(time, FOV, '-k', linewidth = 0.5, label = '$F_{\mathrm{ovS}}$')

ax.plot(time, FOV_mean[0], color = 'r', linestyle = '-', linewidth = 2, label = '1 - 50, 51 - 100, etc.')
ax.scatter(time_average[0], FOV_average[0], color = 'r', marker = 'o', s = 20, zorder = 10, edgecolor = 'firebrick')

ax.plot(time, FOV_mean[24]-0.1, color = 'b', linestyle = '-', linewidth = 2, label = '24 - 74, 75 - 124, etc.')
ax.scatter(time_average[24], FOV_average[24]-0.1, color = 'b', marker = 'D', s = 20, zorder = 10, edgecolor = 'dodgerblue')


ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(1, 2200)
ax.set_ylim(-0.35, 0.35)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = '$F_{\mathrm{ovS}}$')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 1.5, label = '1 - 50, 51 - 100, etc.')
graph_3		= ax.plot([-100, -100], [-100, -100], '-b', linewidth = 1.5, label = '24 - 74, 75 - 124, etc.')

graphs	      = graph_1 + graph_2 + graph_3

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('c) $F_{\mathrm{ovS}}$, cubic splines fit')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax2 	= ax.twinx()

for year_start in range(year_average):
	#Loop over all the cubic splines
	ax2.plot(time, (FOV_mean[year_start]-np.mean(FOV_mean, axis = 0))*1000.0, color = 'gray', linestyle = '-', linewidth = 0.5)

ax.plot(time, FOV, color = 'k', linestyle = '-', linewidth = 0.5)
ax.plot(time, np.mean(FOV_mean, axis = 0), color = 'r', linestyle = '-', linewidth = 2)

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(1, 2200)
ax.set_ylim(-0.35, 0.35)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])

ax.quiver(1732.48, -0.27, 0, 0.4, scale = 5, color = 'r', zorder = 10, width = 0.005)
ax.text(1732.48, -0.29, '$\mathrm{d}_t F_{\mathrm{ovS}} = 0$', verticalalignment='bottom', horizontalalignment='center', color = 'r', fontsize = 9)

#-----------------------------------------------------------------------------------------

ax2.set_ylabel('Freshwater transport difference from $F_{\mathrm{ovS}}$ mean (mSv)')
ax2.set_ylim(-20, 20)

ax2.set_yticks([-20, -10, 0, 10, 20])

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = '$F_{\mathrm{ovS}}$')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 1.5, label = '$F_{\mathrm{ovS}}$ mean')
graph_3		= ax.plot([-100, -100], [-100, -100], '-', color = 'gray', linewidth = 1.5, label = '$F_{\mathrm{ovS}}$ cubic splines')

graphs	      = graph_1 + graph_2 + graph_3

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('d) $F_{\mathrm{ovS}}$, cubic splines')

show()



