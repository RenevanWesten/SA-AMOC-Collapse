#Program determines the AMOC tipping point by using break point regression

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	transport	= fh.variables['Transport'][:]	#AMOC strength (Sv)

	fh.close()

	return time, transport

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

depth_min 	= 0
depth_max	= 1000

#-----------------------------------------------------------------------------------------
time, transport			= ReadinData(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	


tipping_estimate	= []
bins			    = 3
year_bins		    = np.arange(1720, 1800.1, bins)
year_center		    = np.arange(year_bins[0] + bins / 2.0, year_bins[-1], bins)
tipping_point_hist	= np.zeros((len(year_bins)))

for time_1_i in range(1649, 1700):
	#Varying time point 1
	for time_3_j in range(1799, 1850):
		#Varying time point 3

		transport_diff_all	= np.zeros(51)

		for time_2_k in range(1724, 1775):
			#Varying time point 3 (break point)
		
			#Make an empty array for ideal model
			x_break	= np.zeros(time_3_j - time_1_i + 1)

			if transport[time_2_k] > transport[time_1_i]:
                #Positive trend between t2 and t1, do not consider this t2
				transport_diff_all[time_2_k-1724] = 10**9
				continue
            
			#The first part, time 1 to time 2 (break point)
			x_break[:time_2_k-time_1_i+1]	= transport[time_1_i] + np.arange(time_2_k-time_1_i+1) * (transport[time_2_k] - transport[time_1_i]) / (time[time_2_k] - time[time_1_i])

			#The second part, time 2 to time 3
			x_break[time_2_k-time_1_i:]	= transport[time_2_k] + np.arange(time_3_j-time_2_k+1) * (transport[time_3_j] - transport[time_2_k]) / (time[time_3_j] - time[time_2_k])

			#Determine the difference with the real time series
			transport_diff	= transport[time_1_i:time_3_j+1] - x_break
			transport_diff_all[time_2_k-1724]	= np.sum(transport_diff**2.0)

			if time_1_i == 1660 and time_3_j == 1827 and time_2_k == 1753:
				#Minimum for the current settings, save the array
				time_example	= time[time_1_i:time_3_j+1]
				x_break_example	= np.copy(x_break)

		#Save the minimum
		tipping_estimate.append(time[np.argmin(transport_diff_all)+1724])

		index				= np.where(time[np.argmin(transport_diff_all)+1724] >= year_bins)[0][-1]

		tipping_point_hist[index] 	+= 1.0


print('AMOC tipping point (mean):', np.mean(tipping_estimate))
print('AMOC tipping point (5 percentile):', np.percentile(tipping_estimate, 5))
print('AMOC tipping point (10 percentile):', np.percentile(tipping_estimate, 10))
print('AMOC tipping point (25 percentile):', np.percentile(tipping_estimate, 25))
print('AMOC tipping point (75 percentile):', np.percentile(tipping_estimate, 75))
print('AMOC tipping point (90 percentile):', np.percentile(tipping_estimate, 90))
print('AMOC tipping point (95 percentile):', np.percentile(tipping_estimate, 95))

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([1650, 1700], -100, 100, alpha=0.25, edgecolor='orange', facecolor='orange')
ax.fill_between([1725, 1775], -100, 100, alpha=0.25, edgecolor='cyan', facecolor='cyan')
ax.fill_between([1800, 1850], -100, 100, alpha=0.25, edgecolor='red', facecolor='red')

ax.plot(time, transport, '-k')
ax.plot(time_example, x_break_example, '-r', linewidth = 1.5)

ax.plot(time_example[0], x_break_example[0], 'or')
ax.plot(time_example[1753-1660], x_break_example[1753-1660], 'or')
ax.plot(time_example[-1], x_break_example[-1], 'or')

ax.set_xlabel('Model year')
ax.set_ylabel('Volume transport (Sv)')
ax.set_xlim(1600, 1900)
ax.set_ylim(-2, 22)
ax.grid()

ax.text(1675, -1, '$t_1$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 14)
ax.text(1750, -1, '$t_2$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 14)
ax.text(1825, -1, '$t_3$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 14)

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = 'AMOC')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 1.5, label = 'AMOC$_{\mathrm{break}}$')

graphs	      = graph_1 + graph_2

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)


ax.set_title('a) AMOC break regression')

#-----------------------------------------------------------------------------------------

bar_width = 2

fig, ax	= subplots()

ax.bar(year_center, tipping_point_hist[:-1], bar_width, facecolor='white', linewidth = 2, edgecolor = 'royalblue', zorder = 10)
ax.bar(year_center, tipping_point_hist[:-1], bar_width, facecolor='cyan', alpha = 0.25, linewidth = 0, edgecolor = 'royalblue', zorder = 10)

ax.set_xlim(1720, 1780)
ax.set_xlabel('Break point (model year)')
ax.set_ylabel('Frequency')
ax.grid()

ax.set_title('b) AMOC tipping point')
	
show()		

