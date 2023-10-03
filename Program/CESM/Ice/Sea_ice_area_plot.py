#Program plots the Northern and Southern hemispheric sea-ice area time series

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

def SignificantTrend(time, data):
	"""Finds whether trend is significant
	Returns the trend and if it significant (= 1)"""

	#Set time similar to Santer et al. (2000), time array from 1 till N
	#Statistical significance of trends and trend differences in layer-average atmospheric salterature time series
	time		= np.arange(1, len(time) + 1)

	#Determine the detrended time series
	trend, base 	= polyfit(time, data, 1)
	data_res	= data - ((trend * time) + base)

	#Effective sample size, based on the lag-1 correlation
	corr_1		= np.corrcoef(data_res[:-1], data_res[1:])[0, 1]
	N_eff		= int(len(time) * (1.0 - corr_1) / (1.0 + corr_1))

	#Determine the variance of the anomalies
	data_var	= np.sum(data_res**2.0) / (N_eff - 2.0)

	#Determine the standard error
	standard_error	=  np.sqrt(data_var) / np.sqrt(np.sum((time - np.mean(time))**2.0))

	#Determine the Student-T value
	t_value		= trend / standard_error

	#Get the significance levels and the corresponding critical values (two-sided)
	sig_levels 	= np.arange(50, 100, 0.5) / 100.0
	t_crit 		= stats.t.ppf((1.0 + sig_levels) / 2.0, N_eff - 2)

	#Get the indices where the significance is exceeding the critical values
	sig_index	= np.where(fabs(t_value) > t_crit)[0]
	significant	= 0.0

	if len(sig_index) > 0:
		#If there are significance values, take the highest significant level
		significant = sig_levels[sig_index[-1]]

	return trend, np.sqrt(standard_error), significant


#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

fh 		= netcdf.Dataset(directory+'Ice/Sea_ice_area.nc', 'r')

time_ice	= fh.variables['time'][:] 		
ice_area_SH	= fh.variables['AREA_SH'][:] / 10**12 #million km^2 	
ice_area_NH	= fh.variables['AREA_NH'][:] / 10**12 #million km^2 	

fh.close()

fh 		= netcdf.Dataset(directory+'Atmosphere/TEMP_2m_global.nc', 'r')

time		= fh.variables['time'][:] 		
temp_global	= fh.variables['TEMP_global'][:] #2-meter temperature (deg C)
temp_SH		= fh.variables['TEMP_SH'][:] 	#2-meter temperature (deg C)
temp_NH		= fh.variables['TEMP_NH'][:] 	#2-meter temperature (deg C)

fh.close()

#Determine the 2-meter temperature trends and sea-ice area trends
trend_1, error, sig_1 	= SignificantTrend(time[399:1350], temp_global[399:1350])
trend_2, error, sig_2 	= SignificantTrend(time[399:1350], temp_SH[399:1350])
trend_3, error, sig_3 	= SignificantTrend(time[399:1350], temp_NH[399:1350])
trend_4, error, sig_4 	= SignificantTrend(time_ice[4788+2:16200][::12], ice_area_SH[4788+2:16200][::12]+ice_area_NH[4788+2:16200][::12])
trend_5, error, sig_5 	= SignificantTrend(time_ice[4788+2:16200][::12], ice_area_SH[4788+2:16200][::12])
trend_6, error, sig_6 	= SignificantTrend(time_ice[4788+2:16200][::12], ice_area_NH[4788+2:16200][::12])
trend_7, error, sig_7 	= SignificantTrend(time_ice[4788+8:16200][::12], ice_area_SH[4788+8:16200][::12]+ice_area_NH[4788+8:16200][::12])
trend_8, error, sig_8 	= SignificantTrend(time_ice[4788+8:16200][::12], ice_area_SH[4788+8:16200][::12])
trend_9, error, sig_9 	= SignificantTrend(time_ice[4788+8:16200][::12], ice_area_NH[4788+8:16200][::12])

print('Temperature trend over years', int(time[399]), '-',int(time[1349]))
print('Global: '+str(trend_1*100)+', p-value = '+str(1-sig_1))
print('SH: '+str(trend_2*100)+', p-value = '+str(1-sig_2))
print('NH: '+str(trend_3*100)+', p-value = '+str(1-sig_3))
print()

print('Sea-ice area March trend over years', int(time_ice[4788]), '-',int(time_ice[16199]))
print('Global: '+str(trend_4*100)+', p-value = '+str(1-sig_4))
print('SH: '+str(trend_5*100)+', p-value = '+str(1-sig_5))
print('NH: '+str(trend_6*100)+', p-value = '+str(1-sig_6))
print()

print('Sea-ice area September trend over years', int(time_ice[4788]), '-',int(time_ice[16199]))
print('Global: '+str(trend_7*100)+', p-value = '+str(1-sig_7))
print('SH: '+str(trend_8*100)+', p-value = '+str(1-sig_8))
print('NH: '+str(trend_9*100)+', p-value = '+str(1-sig_9))
print()

print('Global')
print('Year', int(time[0]), '-',int(time[49]),':', np.mean(temp_global[0:50]))
print('Year', int(time[2150]), '-',int(time[2199]),':', np.mean(temp_global[2150:2200]), ', Anomaly:', np.mean(temp_global[2150:2200]) - np.mean(temp_global[0:50]))
print()
print('SH')
print('Year', int(time[0]), '-',int(time[49]),':', np.mean(temp_SH[0:50]))
print('Year', int(time[2150]), '-',int(time[2199]),':', np.mean(temp_SH[2150:2200]), ', Anomaly:', np.mean(temp_SH[2150:2200]) - np.mean(temp_SH[0:50]))
print()
print('NH')
print('Year', int(time[0]), '-',int(time[49]),':', np.mean(temp_NH[0:50]))
print('Year', int(time[2150]), '-',int(time[2199]),':', np.mean(temp_NH[2150:2200]), ', Anomaly:', np.mean(temp_NH[2150:2200]) - np.mean(temp_NH[0:50]))

#-----------------------------------------------------------------------------------------

index_march    	= np.arange(2, len(time_ice), 12)
index_september	= np.arange(8, len(time_ice), 12)
time_ice	    = np.arange(len(index_march))+1

year_average 	          = 5
time_average    	    = int(len(time) / year_average)
time_average		    = np.zeros(time_average)
temp_global_average     = np.zeros(len(time_average))
temp_SH_average     	= np.zeros(len(time_average))
temp_NH_average     	= np.zeros(len(time_average))

for time_i in range(len(time_average)):
	#Loop over each time slice
	time_average[time_i]		= np.mean(time[time_i*year_average:(time_i+1)*year_average])
	temp_global_average[time_i] = np.mean(temp_global[time_i*year_average:(time_i+1)*year_average])
	temp_SH_average[time_i] 	= np.mean(temp_SH[time_i*year_average:(time_i+1)*year_average])
	temp_NH_average[time_i] 	= np.mean(temp_NH[time_i*year_average:(time_i+1)*year_average])

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_march	= ax.plot(time_ice, ice_area_NH[index_march], '-c', linewidth = 1.0, label = 'March', zorder = 10)
graph_september	= ax.plot(time_ice, ice_area_NH[index_september], '-b', linewidth = 1.0, label = 'September', zorder = 10)

ax.set_xlabel('Model year')
ax.set_ylabel('Sea-ice area (million km$^2$)')
ax.set_xlim(1, 2200)
ax.set_ylim(0, 30)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])

ax2 = ax.twinx()

graph_global	= ax2.plot(time_average, temp_global_average, linestyle = '-', color = 'k', linewidth = 1.0, label = 'Global')
graph_NH	= ax2.plot(time_average, temp_NH_average, linestyle = '-', color = 'r', linewidth = 1.0, label = 'NH')

ax2.set_ylim(11, 15)
ax2.set_ylabel('2-meter surface temperature ($^{\circ}$C)')
ax2.set_yticks([11, 12, 13, 14, 15])

graphs	      	= graph_march + graph_september

legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc = 'lower left', ncol=1, framealpha = 1.0)

graphs	      	= graph_global + graph_NH

legend_labels 	= [l.get_label() for l in graphs]
legend_2	= ax.legend(graphs, legend_labels, loc = 'lower right', ncol=1, framealpha = 1.0)
ax.add_artist(legend_1)

ax.set_title('c) Northern Hemispheric sea-ice area')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_march	= ax.plot(time_ice, ice_area_SH[index_march], '-c', linewidth = 1.0, label = 'March')
graph_september	= ax.plot(time_ice, ice_area_SH[index_september], '-b', linewidth = 1.0, label = 'September')

ax.set_xlabel('Model year')
ax.set_ylabel('Sea-ice area (million km$^2$)')
ax.set_xlim(1, 2200)
ax.set_ylim(0, 30)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])

ax2 = ax.twinx()

graph_global	= ax2.plot(time_average, temp_global_average, linestyle = '-', color = 'k', linewidth = 1.0, label = 'Global')
graph_SH	= ax2.plot(time_average, temp_SH_average, linestyle = '-', color = 'firebrick', linewidth = 1.0, label = 'SH')

ax2.set_ylim(11, 15)
ax2.set_ylabel('2-meter surface temperature ($^{\circ}$C)')
ax2.set_yticks([11, 12, 13, 14, 15])

graphs	      	= graph_march + graph_september

legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc = 'lower left', ncol=1, framealpha = 1.0)

graphs	      	= graph_global + graph_SH

legend_labels 	= [l.get_label() for l in graphs]
legend_2	= ax.legend(graphs, legend_labels, loc = 'lower right', ncol=1, framealpha = 1.0)
ax.add_artist(legend_1)

ax.set_title('d) Southern Hemispheric sea-ice area')

show()