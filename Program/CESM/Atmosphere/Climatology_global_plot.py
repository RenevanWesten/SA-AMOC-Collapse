#Program plots the global climatology

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
from scipy import stats
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'


def ConverterField(index_break, field):
	"""Shifts field, where it starts at 0E and ends at 360E"""

	new_field	= ma.masked_all(shape(field))
	length_section	= len(field[0]) - index_break

	#Shift the first part
	new_field[:, :length_section] = field[:, index_break:]

	#Shift the last part
	new_field[:, length_section:] = field[:, :index_break] 

	return new_field

def PeriodicBoundaries3D(lon, lat, field, lon_grids = 1):
    """Add periodic zonal boundaries for 2D field"""
    
    #Empty field with additional zonal boundaries
    lon_2                   = np.zeros((len(lat), len(lon[0]) + lon_grids * 2))
    lat_2                   = np.zeros((len(lat), len(lon_2[0])))
    field_2                 = ma.masked_all((len(field), len(lat), len(lon_2[0])))
    
    #Get the left boundary, which is the right boundary of the original field
    lon_2[:, :lon_grids]       	= lon[:, -lon_grids:]
    lat_2[:, :lon_grids]       	= lat[:, -lon_grids:]
    field_2[:, :, :lon_grids]	= field[:, :, -lon_grids:]
    
    #Same for the right boundary
    lon_2[:, -lon_grids:]        	= lon[:, :lon_grids]
    lat_2[:, -lon_grids:]        	= lat[:, :lon_grids]
    field_2[:, :, -lon_grids:]	= field[:, :, :lon_grids]
    
    #And the complete field
    lon_2[:, lon_grids:-lon_grids]          = lon
    lat_2[:, lon_grids:-lon_grids]          = lat
    field_2[:, :, lon_grids:-lon_grids]     = field
    
    return lon_2, lat_2, field_2

def LowCESMPlot(lon, lat, field):
	"""Returns 4 array's to plot on a global projection"""

	#Left of pole
	lon[lon > 180]	= lon[lon > 180] - 360.0

	lon_1		= lon[:, :160]
	lat_1		= lat[:, :160]
	field_1		= field[:, :160]

	#Right of pole
	lon_2		= lon[:, 159:]
	lat_2		= lat[:, 159:]
	field_2		= field[:, 159:]

	lat_3		= ma.masked_where(lon_2 > 0.0, lat_2)
	field_3		= ma.masked_where(lon_2 > 0.0, field_2)
	lon_3		= ma.masked_where(lon_2 > 0.0, lon_2)

	lon_2[lon_2 < -160] 	= lon_2[lon_2 < -160] + 360
	lat_2			= ma.masked_where(lon_2 < 0.0, lat_2)
	field_2			= ma.masked_where(lon_2 < 0.0, field_2)
	lon_2			= ma.masked_where(lon_2 < 0.0, lon_2)

	#To match at 40W
	index_1		= (fabs(lon[40] - 0.0)).argmin()

	lon_4		= ConverterField(index_1, lon)
	lat_4		= ConverterField(index_1, lat)
	field_4		= ConverterField(index_1, field)

	lon_4		= lon_4[:, 280:300]
	lat_4		= lat_4[:, 280:300]
	field_4		= field_4[:, 280:300]

	return lon_1, lat_1, field_1, lon_2, lat_2, field_2, lon_3, lat_3, field_3, lon_4, lat_4, field_4

def Welch(data_1, data_2):
	"""Conducts the Welch t-test"""
	
	#Determine the means
	mean_1	= np.mean(data_1)
	mean_2	= np.mean(data_2)
	
	#Determine the corrected sample standard deviations
	std_1	= np.sqrt(1.0 / (len(data_1) - 1) * np.sum((data_1 - mean_1)**2.0))
	std_2	= np.sqrt(1.0 / (len(data_2) - 1) * np.sum((data_2 - mean_2)**2.0))

	#Determine the Welch t-value
	t_welch	= (mean_1 - mean_2) / np.sqrt((std_1**2.0 / len(data_1)) + (std_2**2.0 / len(data_2)))

	#Determine the degrees of freedome (dof)
	dof	= ((std_1**2.0 / len(data_1)) + (std_2**2.0 / len(data_2)))**2.0 / ((std_1**4.0 / (len(data_1)**2.0 * (len(data_1) - 1))) + (std_2**4.0 / (len(data_2)**2.0 * (len(data_2) - 1))))

	#Get the significance levels and the corresponding critical values (two-sided)
	sig_levels 	= np.arange(50, 100, 0.5) / 100.0
	t_crit		= stats.t.ppf((1.0 + sig_levels) / 2.0, dof)

	#Get the indices where the significance is exceeding the critical values
	sig_index	= np.where(fabs(t_welch) > t_crit)[0]
	significant	= 0.0

	if len(sig_index) > 0:
		#If there are significance values, take the highest significant level
		significant = sig_levels[sig_index[-1]]

	return significant
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

fh 	= netcdf.Dataset(directory+'Atmosphere/Climatology_regions_year_0001-0050.nc', 'r')

temp_1	= fh.variables['TEMP'][:]
prec_1	= fh.variables['PREC'][:] 

fh.close()

fh 	= netcdf.Dataset(directory+'Atmosphere/Climatology_regions_year_2151-2200.nc', 'r')

temp_2	= fh.variables['TEMP'][:]
prec_2	= fh.variables['PREC'][:] 

fh.close()

#Check for significant differences between the climate variables 
diff_sig 	= np.zeros((len(temp_1), len(temp_1[0])))

for region_i in range(len(temp_1)):
	for month_i in range(12):
		#Temperature difference
		temp_sig	= Welch(temp_1[region_i, month_i], temp_2[region_i, month_i])

		#Precipitation difference
		prec_sig	= Welch(prec_1[region_i, month_i], prec_2[region_i, month_i])

		if temp_sig <= 0.95 and prec_sig <= 0.95:
			#Both temperature and precipitation differences are non-significant
			diff_sig[region_i, month_i] = 1

		elif temp_sig <= 0.95:
			#Only temperature difference is non-significant
			diff_sig[region_i, month_i] = 2

		elif prec_sig <= 0.95:
			#Only precipitation difference is non-significant
			diff_sig[region_i, month_i] = 3

#Take the time mean for plotting
temp_1 	= np.mean(temp_1, axis = 2)
temp_2 	= np.mean(temp_2, axis = 2)
prec_1 	= np.mean(prec_1, axis = 2)
prec_2 	= np.mean(prec_2, axis = 2)

month		     = np.arange(1,13)
month_periodic	 = np.arange(0,14)
bar_width 	     = 0.4

#Add periodic boundaries (only for plotting)
temp_1_periodic		= np.zeros((len(temp_1), len(month_periodic)))
temp_2_periodic		= np.zeros((len(temp_2), len(month_periodic)))
temp_1_periodic[:,0]	= temp_1[:, -1]
temp_1_periodic[:,1:13] = temp_1
temp_1_periodic[:,13]	= temp_1[:, 0]
temp_2_periodic[:,0]	= temp_2[:, -1]
temp_2_periodic[:,1:13]	= temp_2
temp_2_periodic[:,13]	= temp_2[:, 0]

#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Robinson()}, figsize = (20, 10))

ax.set_global()
ax.stock_img()
ax.gridlines()
ax.coastlines()

ax.plot([-95, -95, -85, -85, -95], [35, 45, 45, 35, 35], '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax.plot([-70, -70, -60, -60, -70], [-8, 2, 2, -8, -8], '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax.plot([5, 5, 15, 15, 5], [46, 56, 56, 46, 46], '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax.plot([20, 20, 30, 30, 20], [-5, 5, 5, -5, -5], '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax.plot([105, 105, 115, 115, 105], [25, 35, 35, 25, 25], '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax.plot([135, 135, 145, 145, 135], [-30, -20, -20, -30, -30], '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

ax.text(-90, 40, '1', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=16, transform=ccrs.PlateCarree())
ax.text(-65, -3, '2', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=16, transform=ccrs.PlateCarree())
ax.text(10, 51, '3', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=16, transform=ccrs.PlateCarree())
ax.text(25, 0, '4', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=16, transform=ccrs.PlateCarree())
ax.text(110, 30, '5', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=16, transform=ccrs.PlateCarree())
ax.text(140, -25, '6', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=16, transform=ccrs.PlateCarree())
#-----------------------------------------------------------------------------------------
ax1 	= fig.add_axes([0.085, 0.56, 0.15, 0.20])
ax12	= ax1.twinx()

ax1.bar(month - 0.5 * bar_width, prec_1[0], bar_width, color='coral', alpha = 0.5, linewidth = 0.0)
ax1.bar(month + 0.5 * bar_width, prec_2[0], bar_width, color='royalblue', alpha = 0.5, linewidth = 0.0)

ax12.plot(month_periodic, temp_1_periodic[0], '-', linewidth = 2, color = 'r')
ax12.plot(month_periodic, temp_2_periodic[0], '-', linewidth = 2, color = 'b')

ax1.set_xlim(0.5, 12.5)
ax1.grid()
ax1.set_xticks(month)
ax1.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
ax1.set_ylabel('Precipitation (mm)')
ax12.set_ylabel('Temperature ($^{\circ}$C)')

ax1.set_ylim(0, 120)
ax12.set_ylim(-10, 30)


for month_i in range(12):
	#Check for significant differences
	if diff_sig[0, month_i] == 0 or diff_sig[0, month_i] == 3:
		#Temperature is significant different
		ax1.text(month_i / 12 + 0.059, 0.023, 'T', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 7, transform = ax1.transAxes)

	if diff_sig[0, month_i] == 0 or diff_sig[0, month_i] == 2:
		#Precipitation is significant different
		ax1.text(month_i / 12 + 0.025, 0.023, 'P', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 7, transform = ax1.transAxes)

ax1.set_title('1) North America')

#-----------------------------------------------------------------------------------------
ax2 	= fig.add_axes([0.085, 0.27, 0.15, 0.20])
ax22	= ax2.twinx()

ax2.bar(month - 0.5 * bar_width, prec_1[1], bar_width, color='coral', alpha = 0.5, linewidth = 0.0)
ax2.bar(month + 0.5 * bar_width, prec_2[1], bar_width, color='royalblue', alpha = 0.5, linewidth = 0.0)

ax22.plot(month_periodic, temp_1_periodic[1], '-', linewidth = 2, color = 'r')
ax22.plot(month_periodic, temp_2_periodic[1], '-', linewidth = 2, color = 'b')

ax2.set_xlim(0.5, 12.5)
ax2.grid()
ax2.set_xticks(month)
ax2.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
ax2.set_ylabel('Precipitation (mm)')
ax22.set_ylabel('Temperature ($^{\circ}$C)')

ax2.set_ylim(0, 350)
ax22.set_ylim(20, 30)

for month_i in range(12):
	#Check for significant differences
	if diff_sig[1, month_i] == 0 or diff_sig[1, month_i] == 3:
		#Temperature is significant different
		ax2.text(month_i / 12 + 0.059, 0.023, 'T', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 7, transform = ax2.transAxes)

	if diff_sig[1, month_i] == 0 or diff_sig[1, month_i] == 2:
		#Precipitation is significant different
		ax2.text(month_i / 12 + 0.025, 0.023, 'P', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 7, transform = ax2.transAxes)

ax2.set_title('2) Amazon')

#-----------------------------------------------------------------------------------------
ax3 	= fig.add_axes([0.58, 0.7, 0.15, 0.20])
ax32	= ax3.twinx()

ax3.bar(month - 0.5 * bar_width, prec_1[2], bar_width, color='coral', alpha = 0.5, linewidth = 0.0)
ax3.bar(month + 0.5 * bar_width, prec_2[2], bar_width, color='royalblue', alpha = 0.5, linewidth = 0.0)

ax32.plot(month_periodic, temp_1_periodic[2], '-', linewidth = 2, color = 'r')
ax32.plot(month_periodic, temp_2_periodic[2], '-', linewidth = 2, color = 'b')

ax3.set_xlim(0.5, 12.5)
ax3.grid()
ax3.set_xticks(month)
ax3.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
ax3.set_ylabel('Precipitation (mm)')
ax32.set_ylabel('Temperature ($^{\circ}$C)')

ax3.set_ylim(0, 100)
ax32.set_ylim(-20, 20)

for month_i in range(12):
	#Check for significant differences
	if diff_sig[2, month_i] == 0 or diff_sig[2, month_i] == 3:
		#Temperature is significant different
		ax3.text(month_i / 12 + 0.059, 0.023, 'T', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 7, transform = ax3.transAxes)

	if diff_sig[2, month_i] == 0 or diff_sig[2, month_i] == 2:
		#Precipitation is significant different
		ax3.text(month_i / 12 + 0.025, 0.023, 'P', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 7, transform = ax3.transAxes)

ax3.set_title('3) Europe')

#-----------------------------------------------------------------------------------------
ax4 	= fig.add_axes([0.58, 0.25, 0.15, 0.20])
ax42	= ax4.twinx()

ax4.bar(month - 0.5 * bar_width, prec_1[3], bar_width, color='coral', alpha = 0.5, linewidth = 0.0)
ax4.bar(month + 0.5 * bar_width, prec_2[3], bar_width, color='royalblue', alpha = 0.5, linewidth = 0.0)

ax42.plot(month_periodic, temp_1_periodic[3], '-', linewidth = 2, color = 'r')
ax42.plot(month_periodic, temp_2_periodic[3], '-', linewidth = 2, color = 'b')

ax4.set_xlim(0.5, 12.5)
ax4.grid()
ax4.set_xticks(month)
ax4.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
ax4.set_ylabel('Precipitation (mm)')
ax42.set_ylabel('Temperature ($^{\circ}$C)')

ax4.set_ylim(0, 350)
ax42.set_ylim(20, 25)

for month_i in range(12):
	#Check for significant differences
	if diff_sig[3, month_i] == 0 or diff_sig[3, month_i] == 3:
		#Temperature is significant different
		ax4.text(month_i / 12 + 0.059, 0.023, 'T', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 7, transform = ax4.transAxes)

	if diff_sig[3, month_i] == 0 or diff_sig[3, month_i] == 2:
		#Precipitation is significant different
		ax4.text(month_i / 12 + 0.025, 0.023, 'P', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 7, transform = ax4.transAxes)

ax4.set_title('4) Central Africa')

#-----------------------------------------------------------------------------------------
ax5 	= fig.add_axes([0.82, 0.56, 0.15, 0.20])
ax52	= ax5.twinx()

ax5.bar(month - 0.5 * bar_width, prec_1[4], bar_width, color='coral', alpha = 0.5, linewidth = 0.0)
ax5.bar(month + 0.5 * bar_width, prec_2[4], bar_width, color='royalblue', alpha = 0.5, linewidth = 0.0)

ax52.plot(month_periodic, temp_1_periodic[0], '-', linewidth = 2, color = 'r')
ax52.plot(month_periodic, temp_2_periodic[0], '-', linewidth = 2, color = 'b')

ax5.set_xlim(0.5, 12.5)
ax5.grid()
ax5.set_xticks(month)
ax5.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
ax5.set_ylabel('Precipitation (mm)')
ax52.set_ylabel('Temperature ($^{\circ}$C)')

ax5.set_ylim(0, 250)
ax52.set_ylim(-10, 30)

for month_i in range(12):
	#Check for significant differences
	if diff_sig[4, month_i] == 0 or diff_sig[4, month_i] == 3:
		#Temperature is significant different
		ax5.text(month_i / 12 + 0.059, 0.023, 'T', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 7, transform = ax5.transAxes)

	if diff_sig[4, month_i] == 0 or diff_sig[4, month_i] == 2:
		#Precipitation is significant different
		ax5.text(month_i / 12 + 0.025, 0.023, 'P', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 7, transform = ax5.transAxes)

ax5.set_title('5) South-east Asia')

#-----------------------------------------------------------------------------------------
ax6 	= fig.add_axes([0.82, 0.1, 0.15, 0.20])
ax62	= ax6.twinx()

ax6.bar(month - 0.5 * bar_width, prec_1[5], bar_width, color='coral', alpha = 0.5, linewidth = 0.0)
ax6.bar(month + 0.5 * bar_width, prec_2[5], bar_width, color='royalblue', alpha = 0.5, linewidth = 0.0)

ax62.plot(month_periodic, temp_1_periodic[5], '-', linewidth = 2, color = 'r')
ax62.plot(month_periodic, temp_2_periodic[5], '-', linewidth = 2, color = 'b')

ax6.set_xlim(0.5, 12.5)
ax6.grid()
ax6.set_xticks(month)
ax6.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
ax6.set_ylabel('Precipitation (mm)')
ax62.set_ylabel('Temperature ($^{\circ}$C)')

ax6.set_ylim(0, 200)
ax62.set_ylim(10, 30)

for month_i in range(12):
	#Check for significant differences
	if diff_sig[5, month_i] == 0 or diff_sig[5, month_i] == 3:
		#Temperature is significant different
		ax6.text(month_i / 12 + 0.059, 0.023, 'T', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 7, transform = ax6.transAxes)

	if diff_sig[5, month_i] == 0 or diff_sig[5, month_i] == 2:
		#Precipitation is significant different
		ax6.text(month_i / 12 + 0.025, 0.023, 'P', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 7, transform = ax6.transAxes)

ax6.set_title('6) Australia')

#----------------------------------------------------------------------------------------

prec_1_graph 		= mpatches.Patch(facecolor='coral', alpha = 0.5, label='Precipitation (1-50)')
prec_2_graph 		= mpatches.Patch(facecolor='royalblue', alpha = 0.5, label='Precipitation (2151-2200)')
temp_1_graph		= mlines.Line2D([], [], color='r', linewidth = 2.0, label = 'Temperature (1-50)')
temp_2_graph		= mlines.Line2D([], [], color='b', linewidth = 2.0, label = 'Temperature (2151-2200)')

ax.legend(handles=[prec_1_graph, prec_2_graph, temp_1_graph, temp_2_graph], loc=(-0.09, 0), ncol=1, fancybox=True, shadow=False, framealpha = 1.0)

show()