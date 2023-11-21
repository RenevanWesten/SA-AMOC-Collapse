#Program determines the difference in AMOC tipping point and FOV minimum

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

def ReadinDataAMOC(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	transport	= fh.variables['Transport'][:]	#AMOC strength (Sv)

	fh.close()

	return time, transport


def ReadinDataFOV(filename):

	fh = netcdf.Dataset(filename, 'r')

	time	= fh.variables['time'][:]		
	FOV		= fh.variables['F_OV'][:]	#Freshwater transport

	fh.close()

	return time, FOV

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

depth_min 	= 0
depth_max	= 1000

#-----------------------------------------------------------------------------------------
time, transport			= ReadinDataAMOC(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	


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
time, FOV	= ReadinDataFOV(directory+'Ocean/FOV_index_section_34S.nc')
              

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

print()
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

year_diff = np.zeros(1000000)

for diff_i in range(len(year_diff)):
    #Get the difference by drawing random samples
    year_diff[diff_i] = np.random.choice(tipping_estimate) - np.random.choice(FOV_min)
    
print()
print('Difference (mean):', np.mean(year_diff))
print('Difference (5 percentile):', np.percentile(year_diff, 5))
print('Difference (10 percentile):', np.percentile(year_diff, 10))
print('Difference (25 percentile):', np.percentile(year_diff, 25))
print('Difference (75 percentile):', np.percentile(year_diff, 75))
print('Difference (90 percentile):', np.percentile(year_diff, 90))
print('Difference (95 percentile):', np.percentile(year_diff, 95))
