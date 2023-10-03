#Program determines the FOV minimum

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

	time		= fh.variables['time'][:]		
	FOV		    = fh.variables['F_OV'][:]	#Fresh water

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

FOV_mean     = np.zeros((year_average, len(time)))
FOV_der      = np.zeros((year_average, len(time)))

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

#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
time_past, FOV_past = ReadinData(directory+'Ocean/FOV_index_section_34S.nc')

time_past       = time_past[:1696]
FOV_past        = FOV_past[:1696]     

time_average    = int(len(time_past) / year_average)
time_average    = ma.masked_all((year_average, time_average))
FOV_average     = ma.masked_all((year_average, len(time_average[0])))

for year_start in range(year_average):
    for time_i in range(len(time_average[0])):
        #Loop over each time slice
        if len(time_past[year_start+time_i*year_average:year_start+(time_i+1)*year_average]) != year_average:
            #No complete period
            continue
            
        time_average[year_start, time_i] = np.mean(time_past[year_start+time_i*year_average:year_start+(time_i+1)*year_average])    
        FOV_average[year_start, time_i]  = np.mean(FOV_past[year_start+time_i*year_average:year_start+(time_i+1)*year_average])    

#-----------------------------------------------------------------------------------------	

time_2        = np.arange(1500, time_past[-1])
FOV_der_2     = np.zeros((year_average, len(time_2)))

for year_start in range(year_average):
    #Get the different knots
    time_knot       = time_average[year_start]
    FOV_knot        = FOV_average[year_start]

    #Get non-masked data
    mask_index      = np.where(time_knot.mask == False)[0]
    time_knot       = time_knot[mask_index]
    FOV_knot        = FOV_knot[mask_index]

    #Fit the cubic splines
    cs                    = CubicSpline(time_knot, FOV_knot, bc_type = 'natural')
    FOV_der_2[year_start] = cs(time_2, 1)
    
trend, base 	= polyfit(time_2, np.mean(FOV_der_2, axis = 0), 1)

print(-base / trend)

#-----------------------------------------------------------------------------------------

fig, ax = subplots()

ax.fill_between([0, 1600], -10, 10, alpha=0.25, edgecolor='orange', facecolor='orange')
ax.fill_between([1600, 1695], -10, 10, alpha=0.25, edgecolor='cyan', facecolor='cyan')
ax.fill_between([1695, 1900], -10, 10, alpha=0.25, edgecolor='red', facecolor='red')

for year_start in range(year_average):
    ax.plot(time_2, FOV_der_2[year_start] * 1000, '-', color = 'gray', linewidth = 0.5)

ax.plot(time_2, np.mean(FOV_der_2, axis = 0) * 1000, '-r', linewidth = 2.0)
ax.plot(time[1694:], np.mean(FOV_der[:, 1694:], axis = 0) * 1000, ':', color = 'firebrick', linewidth = 2.0)
ax.plot(time_2, (time_2 * trend + base) * 1000, '-b', linewidth = 2.0)
ax.plot(time[1694:], (time[1694:] * trend + base) * 1000, ':', color = 'dodgerblue', linewidth = 2.0)

ax.set_xlabel('Model year')
ax.set_ylabel(r'Freshwater transport derivative (mSv yr$^{-1}$)')
ax.set_xlim(1500, 1800)
ax.set_ylim(-0.5, 0.5)
ax.grid()

ax.scatter(1745.4172686852426, 0, s = 50, color = 'b', marker = 'D', edgecolor = 'dodgerblue', zorder = 10)

graph_1 = ax.plot([0,0], [-10,-10], '-', color = 'gray', linewidth = 2, label = '$\mathrm{d}_t F_{\mathrm{ovS}}$ (1 - 1695)')
graph_2 = ax.plot([0,0], [-10,-10], '-', color = 'r', linewidth = 2, label = '$\mathrm{d}_t F_{\mathrm{ovS}}$ mean (1 - 1695)')
graph_3 = ax.plot([0,0], [-10,-10], ':', color = 'firebrick', linewidth = 2, label = '$\mathrm{d}_t F_{\mathrm{ovS}}$ mean (1 - 2200)')
graph_4 = ax.plot([0,0], [-10,-10], '-', color = 'b', linewidth = 2, label = 'Linear fit (1500 - 1695)')
graph_5 = ax.plot([0,0], [-10,-10], ':', color = 'dodgerblue', linewidth = 2, label = 'Linear fit (extrapolated)')


graphs	      = graph_1 + graph_2 + graph_3 + graph_4 + graph_5

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.text(0.17, 0.035, 'Available data', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 10, transform = ax.transAxes)
ax.text(0.49, 0.035, 'Future data', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 10, transform = ax.transAxes)
ax.text(0.82, 0.035, 'Extrapolation', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 10, transform = ax.transAxes)

ax.set_title('a) $F_{\mathrm{ovS}}$ minimum estimate procedure')      

#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

year_averages_all   = [40, 45, 50]
year                = np.arange(1600, 1733)
FOV_min             = np.zeros((len(year_averages_all), len(year)))

for average_i in range(len(year_averages_all)):
    year_average    = int(year_averages_all[average_i])
    for year_i in range(len(year)):
        #Read in the time series of the PI control and forcing (historical + RCP8.5)
        time_past, FOV_past	= ReadinData(directory+'Ocean/FOV_index_section_34S.nc')
        
        year_end     = int(year[year_i])
        time_past    = time_past[:year_end]
        FOV_past     = FOV_past[:year_end]
        
        time_average    = int(len(time) / year_average)
        time_average    = ma.masked_all((year_average, time_average))
        FOV_average     = ma.masked_all((year_average, len(time_average[0])))
        
        for year_start in range(year_average):
            for time_i in range(len(time_average[0])):
                #Loop over each time slice
                if len(time_past[year_start+time_i*year_average:year_start+(time_i+1)*year_average]) != year_average:
                    #No complete period
                    continue
                    
            
                time_average[year_start, time_i] = np.mean(time_past[year_start+time_i*year_average:year_start+(time_i+1)*year_average])    
                FOV_average[year_start, time_i]  = np.mean(FOV_past[year_start+time_i*year_average:year_start+(time_i+1)*year_average])    
        
        
        time_2        = np.arange(1500, time_past[-1])
        FOV_der_2     = np.zeros((year_average, len(time_2)))
        
        for year_start in range(year_average):
            #Get the different knots
            time_knot       = time_average[year_start]
            FOV_knot        = FOV_average[year_start]
        
            #Get non-masked data
            mask_index      = np.where(time_knot.mask == False)[0]
            time_knot       = time_knot[mask_index]
            FOV_knot        = FOV_knot[mask_index]
        
            #Fit the cubic splines
            cs                    = CubicSpline(time_knot, FOV_knot, bc_type = 'natural')
            FOV_der_2[year_start] = cs(time_2, 1)
            
        trend, base 	= polyfit(time_2, np.mean(FOV_der_2, axis = 0), 1)
        FOV_min[average_i, year_i] = -base / trend

#-----------------------------------------------------------------------------------------

fig, ax = subplots()

graph_1 = ax.plot(year, FOV_min[0], '-k', linewidth = 2.0, label = '40-year averages')
graph_2 = ax.plot(year, FOV_min[1], '-r', linewidth = 2.0, label = '45-year averages')
graph_3 = ax.plot(year, FOV_min[2], '-b', linewidth = 2.0, label = '50-year averages')
graph_4 = ax.plot([1400, 1900], [1732.48, 1732.48], ':', color = 'firebrick', linewidth = 2.0, label = '$\mathrm{d}_t F_{\mathrm{ovS}} = 0$')
graph_4 = ax.plot([1732.48, 1732.48], [1400, 2100], ':', color = 'firebrick', linewidth = 2.0, label = '$\mathrm{d}_t F_{\mathrm{ovS}} = 0$')

ax.scatter(1695, 1745.4172686852426, s = 50, color = 'b', marker = 'D', edgecolor = 'dodgerblue', zorder = 10)

ax.set_xlabel('Model years included in analysis')
ax.set_ylabel('Model year estimate to $F_{\mathrm{ovS}}$ minimum')
ax.set_xlim(1500, 1800)
ax.set_ylim(1600, 2000)
ax.grid()

graphs	      = graph_1 + graph_2 + graph_3 + graph_4

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)


ax.set_title('b) $F_{\mathrm{ovS}}$ minimum estimate')      

show()


