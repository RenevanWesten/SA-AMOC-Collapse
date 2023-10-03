#Program determines the 2-meter surface temperature trend

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
from scipy import stats

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'


def ReadinData(filename):
	"""Read in the data"""
	fh = netcdf.Dataset(filename, 'r')

	year 		= fh.variables['year'][:]		#Model year
	lon 		= fh.variables['lon'][:]		#Longitude
	lat 		= fh.variables['lat'][:]		#Latitude 
	temp		= fh.variables['TEMP_2m'][:]    #Temperature (deg C)

	fh.close()
        
	return year, lon, lat, temp

def SignificantTrend(time, data):
	"""Finds whether trend is significant
	Returns the trend and if it significant (= 1)"""

	#Set time similar to Santer et al. (2000), time array from 1 till N
	#Statistical significance of trends and trend differences in layer-average atmospheric temperature time series
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

year_start	= 1750
year_end	= 1850

files = glob.glob(directory+'Data/TEMP_2m/CESM_year_*.nc')
files.sort()

#-----------------------------------------------------------------------------------------

#Determine the section length per depth layer
year, lon, lat, temp	= ReadinData(files[0])

time_all        = np.zeros(len(files) * len(year))
temp_all		= ma.masked_all((len(time_all), len(lat), len(lon)))


for file_i in range(len(files)):
    #Now determine for AMOC strength
    print(file_i)
	    
    year, lon, lat, temp = ReadinData(files[file_i])

    for year_i in range(len(year)):
        #Determine the meridional transport
        time_all[file_i*len(year)+year_i] = year[year_i]
        temp_all[file_i*len(year)+year_i] = temp[year_i]

time_start	= (np.abs(time_all - year_start)).argmin()
time_end	= (np.abs(time_all - (year_end+1))).argmin()

time_all	= time_all[time_start:time_end]
temp_all	= temp_all[time_start:time_end]

#-----------------------------------------------------------------------------------------
#Define empty fields
temp_trend		    = ma.masked_all((len(lat), len(lon)))
temp_trend_sig		= ma.masked_all((len(lat), len(lon)))

#-----------------------------------------------------------------------------------------

for lat_i in range(len(lat)):
	print(lat_i, len(lat))
	for lon_i in range(len(lon)):

		trend, error, sig = SignificantTrend(time_all, temp_all[:, lat_i, lon_i])

		#Determine trend per century
		temp_trend[lat_i, lon_i]	= trend * 100.0

		#Save the significance
		temp_trend_sig[lat_i, lon_i]	= sig

#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'Atmosphere/TEMP_2m_trend_year_'+str(year_start)+'-'+str(year_end)+'.nc', 'w')

fh.createDimension('lat', len(lat))
fh.createDimension('lon', len(lon))

fh.createVariable('lat', float, ('lat'), zlib=True)
fh.createVariable('lon', float, ('lon'), zlib=True)
fh.createVariable('TEMP_trend', float, ('lat', 'lon'), zlib=True)
fh.createVariable('TEMP_trend_sig', float, ('lat', 'lon'), zlib=True)

fh.variables['lon'].longname 		= 'Array of longitudes'
fh.variables['lat'].longname 		= 'Array of latitudes'
fh.variables['TEMP_trend'].longname 	= '2-meter surface temperature trend'
fh.variables['TEMP_trend_sig'].longname = 'Level of significance'

fh.variables['lon'].units 		= 'Degrees E'
fh.variables['lat'].units 		= 'Degrees N'
fh.variables['TEMP_trend'].units 	= 'deg C per century'

#Writing data to correct variable	
fh.variables['lon'][:] 			= lon
fh.variables['lat'][:] 			= lat
fh.variables['TEMP_trend'][:] 		= temp_trend
fh.variables['TEMP_trend_sig'][:] 	= temp_trend_sig

fh.close()
