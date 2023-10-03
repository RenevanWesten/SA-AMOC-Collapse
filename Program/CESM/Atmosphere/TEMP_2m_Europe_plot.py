#Program plots the 2-meter surface temperature for various European cities

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

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'


def YearlyConverter(time, data, month_start = 1, month_end = 12, yearly_sum = False):
	"""Determines yearly averaged, over different months of choice,
	default is set to January - December"""

	#Take twice the amount of years for the month day
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31., 31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= month_days[month_start - 1:month_end]
	month_days	= month_days / np.sum(month_days)
	
	if month_end <= 12:
		#Normal average over a single year, for example, February 100 - December 100
		time_year		= np.zeros(int(len(time) / 12))

	else:
		#If you take the average, for example, over November 100 - May 101
		#Take year 101 as the average over this period
		#There is one year less compared to the period analysed
		time_year		= np.zeros(int(len(time) / 12) - 1)

	#-----------------------------------------------------------------------------------------
	data_year	= ma.masked_all(len(time_year))

	for year_i in range(len(time_year)):
		#Determine the SSH over the selected months

		#The year is defined as the current year
		year			= int(time[year_i * 12])

		if month_end	>= 13:
			#If average is taken over, for example, November 100 - May 101, the year is defined as 101
			year = year + 1

		time_year[year_i] 	= year

		#Determine the time mean over the months of choice
		data_year[year_i]		= np.sum(data[year_i * 12 + month_start - 1: year_i * 12 + month_end] * month_days, axis = 0)

		if yearly_sum:
			#Take the yearly sum over the months of choice
			data_year[year_i]	= np.sum(data[year_i * 12 + month_start - 1: year_i * 12 + month_end], axis = 0)

	return time_year, data_year




#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

depth_min 	= 0
depth_max	= 1000

fh = netcdf.Dataset(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc', 'r')

time_AMOC	= fh.variables['time'][:]		
transport	= fh.variables['Transport'][:]	#AMOC strength (Sv)

fh.close()

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Atmosphere/TEMP_2m_Europe_locations.nc', 'r')

time		= fh.variables['time'][:]		
temp		= fh.variables['TEMP'][:]	#Temperature for 5 locations

fh.close()

#-----------------------------------------------------------------------------------------

time_year, temp_1	= YearlyConverter(time, temp[:, 0])
time_year, temp_2	= YearlyConverter(time, temp[:, 1])
time_year, temp_3	= YearlyConverter(time, temp[:, 2])
time_year, temp_4	= YearlyConverter(time, temp[:, 3])
time_year, temp_5	= YearlyConverter(time, temp[:, 4])

#Determine the trends for each month
temp_trend	= np.zeros((5, 14))

year_start	= 1750
year_end	= 1850

time_start	= (np.abs(time - year_start)).argmin()
time_end	= (np.abs(time - (year_end+1))).argmin()

time		= time[time_start:time_end]
temp		= temp[time_start:time_end]

for month_i in range(14):
	#Loop over each month
	if month_i == 0:
		#December
		time_index	= np.arange(11, len(time), 12)
	elif month_i == 13:
		#January
		time_index	= np.arange(0, len(time), 12)
	else:
		#Standard months
		time_index	= np.arange(month_i-1, len(time), 12)

	for loc_i in range(5):
		#Determine the trend per month for each location
		trend, base 			= polyfit(np.arange(len(time_index)), temp[time_index, loc_i], 1)
		temp_trend[loc_i, month_i]	= trend * 100.

#-----------------------------------------------------------------------------------------
fig, ax	= subplots()

ax.fill_between([1750, 1850], -30, 10, alpha=0.25, edgecolor='orange', facecolor='orange')

graph_london	= plot(time_year, temp_1 - temp_1[0], linestyle = '-', color = 'k', linewidth = 1.5, label = 'London')
graph_madrid	= plot(time_year, temp_2 - temp_2[0], linestyle = '-', color = 'r', linewidth = 1.5, label = 'Madrid')
graph_vienna	= plot(time_year, temp_3 - temp_3[0], linestyle = '-', color = 'b', linewidth = 1.5, label = 'Vienna')
graph_bergen	= plot(time_year, temp_4 - temp_4[0], linestyle = '-', color = 'c', linewidth = 1.5, label = 'Bergen')
graph_reykjavik	= plot(time_year, temp_5 - temp_5[0], linestyle = '-', color = 'firebrick', linewidth = 1.5, label = 'Reykjavik')

ax.set_xlabel('Model year')
ax.set_ylabel('Temperature difference ($^{\circ}$C)')
ax.set_xlim(1600, 2000)
ax.set_ylim(-20, 5)
ax.grid()

ax2 	= ax.twinx()

graph_AMOC	= ax2.plot(time_AMOC, transport, linestyle = '-', color = 'gray', linewidth = 1.5, label = 'AMOC')

ax2.set_ylim(-2, 22)
ax2.set_ylabel('Volume transport (Sv)')

#-----------------------------------------------------------------------------------------

graphs	      = graph_london + graph_madrid + graph_vienna + graph_bergen + graph_reykjavik + graph_AMOC

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('c) Yearly temperature and AMOC strength')

#-----------------------------------------------------------------------------------------

month		= np.arange(0,14)

fig, ax	= subplots()

graph_london	= plot(month, temp_trend[0], linestyle = '-', color = 'k', linewidth = 2, label = 'London')
graph_madrid	= plot(month, temp_trend[1], linestyle = '-', color = 'r', linewidth = 2, label = 'Madrid')
graph_vienna	= plot(month, temp_trend[2], linestyle = '-', color = 'b', linewidth = 2, label = 'Vienna')
graph_bergen	= plot(month, temp_trend[3], linestyle = '-', color = 'c', linewidth = 2, label = 'Bergen')
graph_reykjavik	= plot(month, temp_trend[4], linestyle = '-', color = 'firebrick', linewidth = 2, label = 'Reykjavik')

ax.set_xlim(0.5, 12.5)
ax.set_ylim(-40, 5)
ax.grid()
ax.set_xticks(np.arange(1,13))
ax.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
ax.set_ylabel('Temperature trend ($^{\circ}$C per century)')

#-----------------------------------------------------------------------------------------
graphs	      = graph_london + graph_madrid + graph_vienna + graph_bergen + graph_reykjavik

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='lower center', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('d) Monthly temperature trend (1750 - 1850)')

show()
