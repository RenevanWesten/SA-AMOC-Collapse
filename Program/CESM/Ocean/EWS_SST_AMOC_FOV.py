#Program determines the EWS of the subpolar region, AMOC strength and FOV

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

def EWS(time_1, time_2, time_all, data_all, window = 70):
	"""Determines the variance lag-1 auto-correlation following Boers (2021)"""

	time_index_1	= (np.abs(time_all - time_1)).argmin()
	time_index_2	= (np.abs(time_all - time_2)).argmin()+1
	time, data	= time_all[time_index_1:time_index_2], data_all[time_index_1:time_index_2]

	#-----------------------------------------------------------------------------------------	

	auto_corr	= ma.masked_all(len(time))
	var		= ma.masked_all(len(time))

	trend, base	= polyfit(time, data, 1)
	data		= data - ((trend * time) + base)

	for time_i in range(len(time) - window + 1):
		#Loop over each time window
		data_window	= data[time_i:time_i+window]

		#Now determine the auto-correlation and variance
		auto_corr[time_i + int(window / 2)]	= np.corrcoef(data_window[:-1], data_window[1:])[0,1]
		var[time_i + int(window / 2)]		= np.var(data_window)

	return time, trend, base, data, auto_corr, var

def YearlyConverter(time, data, month_start = 1, month_end = 12):
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
		time_year		= np.zeros(int(len(time) / 12 - 1))
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

	return time_year, data_year


def Fourrier_surrogates(data, num_surr):
	"""Takes the fourier transorm of the time series and re-shuffles the statistics"""
	data_fourier        = np.fft.rfft(data)
	random_phases       = np.exp(np.random.uniform(0, 2 * np.pi, (num_surr, len(data) // 2 + 1)) * 1.0j)
	data_fourier_new    = data_fourier * random_phases
	data_new            = np.real(np.fft.irfft(data_fourier_new))

	return data_new

def Kendall_tau_test(data, num_surr, time_data = False):
	"""Conducts the Kendall-Fourier test (adapted from Boers 2021)"""

	#Remove masked elements (if present)
	mask_index  = np.where(data.mask == False)[0]
	data	    = data[mask_index]

	if type(time_data) == type(False):
		#No time array is provided, use dummy variable
		trend, base = stats.linregress(np.arange(len(data)),data)[0], stats.linregress(np.arange(len(data)),data)[1]
	else:
		#Time array is provided, use these to fit the trends
		time_data	= time_data[mask_index]
		trend, base 	= stats.linregress(time_data,data)[0], stats.linregress(time_data,data)[1]	

	#Make time mean zero and create surrogate data
	data_0      = data - np.mean(data)
	data_surr   = Fourrier_surrogates(data_0, num_surr)

	#Determine the linear regression of the surrogate time series
	stat_surr   = np.zeros(num_surr)

	for surr_i in range(num_surr):
		#Determine the linear regression
		stat_surr[surr_i] = stats.linregress(np.arange(len(data)), data_surr[surr_i])[0]

	#Determine the significant level w.r.t. trend
	p = 1 - stats.percentileofscore(stat_surr, trend) / 100.

	return p, trend, base


#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

window		= 70 	#Sliding window (years)

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Ocean/SST_subpolar.nc', 'r')

time		    = fh.variables['time'][:]     	  	
temp_subpolar	= fh.variables['SST_sub'][:] 		
temp_global	    = fh.variables['SST_global'][:]		

fh.close()

time_year, temp_subpolar_all	= YearlyConverter(time, temp_subpolar, 11, 17)
time_all, temp_global		    = YearlyConverter(time, temp_global, 11, 17)

temp_subpolar_all		= temp_subpolar_all - temp_global
temp_subpolar_all		= temp_subpolar_all - np.mean(temp_subpolar_all[:50])

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Ocean/AMOC_transport_depth_0-1000m.nc', 'r')

time_AMOC	= fh.variables['time'][:]		
transport	= fh.variables['Transport'][:]	#AMOC strength (Sv)

fh.close()
#-----------------------------------------------------------------------------------------


fh = netcdf.Dataset(directory+'Ocean/FOV_index_section_34S.nc', 'r')

time_FOV	= fh.variables['time'][:]		
FOV_34S		= fh.variables['F_OV'][:]	#FOV at 34S (Sv)

fh.close()

#-----------------------------------------------------------------------------------------

time_1, trend_1, base_1, data_1, auto_corr_1, var_1	= EWS(600, 900, time_all, temp_subpolar_all, window)
time_2, trend_2, base_2, data_2, auto_corr_2, var_2	= EWS(1000, 1300, time_all, temp_subpolar_all, window)
time_3, trend_3, base_3, data_3, auto_corr_3, var_3	= EWS(1400, 1700, time_all, temp_subpolar_all, window)
time_4, trend_4, base_4, data_4, auto_corr_4, var_4	= EWS(600, 900, time_AMOC, transport, window)
time_5, trend_5, base_5, data_5, auto_corr_5, var_5	= EWS(1000, 1300, time_AMOC, transport, window)
time_6, trend_6, base_6, data_6, auto_corr_6, var_6	= EWS(1400, 1700, time_AMOC, transport, window)
time_7, trend_7, base_7, data_7, auto_corr_7, var_7	= EWS(600, 900, time_FOV, FOV_34S, window)
time_8, trend_8, base_8, data_8, auto_corr_8, var_8	= EWS(1000, 1300, time_FOV, FOV_34S, window)
time_9, trend_9, base_9, data_9, auto_corr_9, var_9	= EWS(1400, 1700, time_FOV, FOV_34S, window)

#Determine the shuffled Fourier Kendall test (following Boers 2021)
p_auto_corr_1, trend_auto_corr_1, base_auto_corr_1	= Kendall_tau_test(auto_corr_1, 100000, time_1)
p_auto_corr_2, trend_auto_corr_2, base_auto_corr_2	= Kendall_tau_test(auto_corr_2, 100000, time_2)
p_auto_corr_3, trend_auto_corr_3, base_auto_corr_3	= Kendall_tau_test(auto_corr_3, 100000, time_3)
p_auto_corr_4, trend_auto_corr_4, base_auto_corr_4	= Kendall_tau_test(auto_corr_4, 100000, time_4)
p_auto_corr_5, trend_auto_corr_5, base_auto_corr_5	= Kendall_tau_test(auto_corr_5, 100000, time_5)
p_auto_corr_6, trend_auto_corr_6, base_auto_corr_6	= Kendall_tau_test(auto_corr_6, 100000, time_6)
p_auto_corr_7, trend_auto_corr_7, base_auto_corr_7	= Kendall_tau_test(auto_corr_7, 100000, time_7)
p_auto_corr_8, trend_auto_corr_8, base_auto_corr_8	= Kendall_tau_test(auto_corr_8, 100000, time_8)
p_auto_corr_9, trend_auto_corr_9, base_auto_corr_9	= Kendall_tau_test(auto_corr_9, 100000, time_9)
p_var_1, trend_var_1, base_var_1			= Kendall_tau_test(var_1, 100000, time_1)
p_var_2, trend_var_2, base_var_2			= Kendall_tau_test(var_2, 100000, time_2)
p_var_3, trend_var_3, base_var_3			= Kendall_tau_test(var_3, 100000, time_3)
p_var_4, trend_var_4, base_var_4			= Kendall_tau_test(var_4, 100000, time_4)
p_var_5, trend_var_5, base_var_5			= Kendall_tau_test(var_5, 100000, time_5)
p_var_6, trend_var_6, base_var_6			= Kendall_tau_test(var_6, 100000, time_6)
p_var_7, trend_var_7, base_var_7			= Kendall_tau_test(var_7, 100000, time_7)
p_var_8, trend_var_8, base_var_8			= Kendall_tau_test(var_8, 100000, time_8)
p_var_9, trend_var_9, base_var_9			= Kendall_tau_test(var_9, 100000, time_9)

#-----------------------------------------------------------------------------------------	
trend, base 	= polyfit(time_all, temp_subpolar_all, 1)
	
fig, ax	= subplots()

ax.fill_between([time_1[0], time_1[-1]], -100, 100, alpha=0.25, edgecolor='orange', facecolor='orange')
ax.fill_between([time_2[0], time_2[-1]], -100, 100, alpha=0.25, edgecolor='cyan', facecolor='cyan')
ax.fill_between([time_3[0], time_3[-1]], -100, 100, alpha=0.25, edgecolor='red', facecolor='red')

graph_1	= ax.plot(time_all, temp_subpolar_all, '-b', linewidth = 0.5, label = 'SST subpolar region')

ax.plot(time_1, ((trend_1 * time_1) + base_1), '-r', linewidth = 2)
ax.plot(time_2, ((trend_2 * time_2) + base_2), '-r', linewidth = 2)
ax.plot(time_3, ((trend_3 * time_3) + base_3), '-r', linewidth = 2)

ax.set_xlabel('Model year')
ax.set_ylabel('Temperature differene ($^{\circ}$C)', color = 'b')
ax.set_xlim(550, 1850)
ax.set_ylim(-6, 6)
ax.grid()

for tl in ax.get_yticklabels():
    tl.set_color('b')

ax2 	= ax.twinx()
graph_2	= ax2.plot(time_AMOC, transport, '-k', linewidth = 0.5, label = 'AMOC strength')
ax2.plot(time_4, ((trend_4 * time_4) + base_4), '-r', linewidth = 2)
ax2.plot(time_5, ((trend_5 * time_5) + base_5), '-r', linewidth = 2)
ax2.plot(time_6, ((trend_6 * time_6) + base_6), '-r', linewidth = 2)
ax2.set_ylabel('Volume transport (Sv)')
ax2.set_ylim(-2, 22)


graph_1	= ax.plot([-100, -100], [-20, -20], '-b', linewidth = 1.0, label = 'SST subpolar region')
graph_2	= ax.plot([-100, -100], [-20, -20], '-k', linewidth = 1.0, label = 'AMOC strength')
graph_3	= ax.plot([-100, -100], [-20, -20], '-r', linewidth = 1.0, label = 'Linear fit')

graphs	      = graph_1 + graph_2 + graph_3

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.text(0.155, 0.035, '$F_{\mathrm{ovS}} > 0$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.465, 0.035, r'$F_{\mathrm{ovS}} \approx 0$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.77, 0.035, r'$F_{\mathrm{ovS}} < 0$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)

ax.set_title('a) SST subpolar region')

#-----------------------------------------------------------------------------------------
ax2 	= fig.add_axes([0.59, 0.60, 0.35, 0.35], projection = ccrs.NearsidePerspective(-30, 50, 10000000))

ax2.coastlines(resolution='110m')
ax2.gridlines()
ax2.add_feature(cfeature.LAND, zorder=10)
ax2.set_global()


lon_domain_sub	= [300, 306, 312, 313, 314, 318, 319, 320, 322, 323, 324, 326, 327, 328, 329, 335, 337, 343, 347, 347, 345, 344, 343, 342, 339, 339, 338, 337, 337, 336, 335, 333, 332, 322, 321, 319, 318, 316, 314, 313, 313, 312, 311, 308, 308, 307, 306, 303, 302, 300, 299, 297, 296, 296, 297, 298, 298, 299, 300]
lat_domain_sub	= [66, 66, 60, 60, 59, 59, 60, 60, 62, 62, 63, 63, 64, 64, 65, 65, 63, 63, 59, 53, 51, 51, 50, 50, 47, 46, 46, 45, 44, 44, 43, 43, 42, 42, 43, 43, 44, 44, 46, 46, 47, 48, 48, 51, 52, 53, 53, 56, 56, 58, 58, 60, 60, 63, 63, 64, 65, 65, 66]

ax2.fill(lon_domain_sub, lat_domain_sub, facecolor='lightskyblue', edgecolor='dodgerblue', linewidth=2, transform=ccrs.PlateCarree())

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([time_1[0], time_1[-1]], -100, 100, alpha=0.25, edgecolor='orange', facecolor='orange')
ax.fill_between([time_2[0], time_2[-1]], -100, 100, alpha=0.25, edgecolor='cyan', facecolor='cyan')
ax.fill_between([time_3[0], time_3[-1]], -100, 100, alpha=0.25, edgecolor='red', facecolor='red')

ax.plot(time_1, trend_var_1 * time_1 + base_var_1, '--', color = 'gray', linewidth = 1)
ax.plot(time_2, trend_var_2 * time_2 + base_var_2, '--', color = 'gray', linewidth = 1)
ax.plot(time_3, trend_var_3 * time_3 + base_var_3, '--', color = 'gray', linewidth = 1)

graph_1	= ax.plot(time_1, var_1, '-k', linewidth = 1, label = 'Variance')
graph_1	= ax.plot(time_2, var_2, '-k', linewidth = 1, label = 'Variance')
graph_1	= ax.plot(time_3, var_3, '-k', linewidth = 1, label = 'Variance')

ax.set_xlabel('Model year')
ax.set_ylabel('Variance ($^{\circ}$C$^2$)')
ax.set_xlim(550, 1850)
ax.set_ylim(0, 0.2)
ax.set_yticks([0, 0.05, 0.10, 0.15, 0.20])
ax.grid()

ax2 = ax.twinx()

ax2.plot(time_1, trend_auto_corr_1 * time_1 + base_auto_corr_1, '--', color = 'dodgerblue', linewidth = 1)
ax2.plot(time_2, trend_auto_corr_2 * time_2 + base_auto_corr_2, '--', color = 'dodgerblue', linewidth = 1)
ax2.plot(time_3, trend_auto_corr_3 * time_3 + base_auto_corr_3, '--', color = 'dodgerblue', linewidth = 1)

graph_2	= ax2.plot(time_1, auto_corr_1, '-b', linewidth = 1, label = 'Auto-correlation')
graph_2	= ax2.plot(time_2, auto_corr_2, '-b', linewidth = 1, label = 'Auto-correlation')
graph_3	= ax2.plot(time_3, auto_corr_3, '-b', linewidth = 1, label = 'Auto-correlation')

ax2.set_ylabel('Lag-1 auto-correlation', color = 'b')
ax2.set_ylim(0, 1)

for tl in ax2.get_yticklabels():
    tl.set_color('b')

graphs	      = graph_1 + graph_2

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('b) SST subpolar region, variance and auto-correlation')

ax.text(0.155, 0.035, '$F_{\mathrm{ovS}} > 0$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.465, 0.035, r'$F_{\mathrm{ovS}} \approx 0$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.77, 0.035, r'$F_{\mathrm{ovS}} < 0$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)

#Significance of trends (period 1)
ax.text(0.155, 0.2, '$p > 0.1$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.155, 0.15, '$p > 0.1$', verticalalignment='center', horizontalalignment='center', color = 'b', fontsize = 12, transform = ax.transAxes)

#Significance of trends (period 2)
ax.text(0.465, 0.2, '$p > 0.1$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.465, 0.15, '$p > 0.1$', verticalalignment='center', horizontalalignment='center', color = 'b', fontsize = 12, transform = ax.transAxes)

#Significance of trends (period 3)
ax.text(0.77, 0.2, '$p > 0.1$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.77, 0.15, '$p > 0.1$', verticalalignment='center', horizontalalignment='center', color = 'b', fontsize = 12, transform = ax.transAxes)

print(p_var_1, p_var_2, p_var_3)
print(p_auto_corr_1, p_auto_corr_2, p_auto_corr_3)

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([time_4[0], time_4[-1]], -100, 100, alpha=0.25, edgecolor='orange', facecolor='orange')
ax.fill_between([time_5[0], time_5[-1]], -100, 100, alpha=0.25, edgecolor='cyan', facecolor='cyan')
ax.fill_between([time_6[0], time_6[-1]], -100, 100, alpha=0.25, edgecolor='red', facecolor='red')

ax.plot(time_4, trend_var_4 * time_4 + base_var_4, '--', color = 'gray', linewidth = 1)
ax.plot(time_5, trend_var_5 * time_5 + base_var_5, '--', color = 'gray', linewidth = 1)
ax.plot(time_6, trend_var_6 * time_6 + base_var_6, '--', color = 'gray', linewidth = 1)

graph_1	= ax.plot(time_4, var_4, '-k', linewidth = 1, label = 'Variance')
graph_1	= ax.plot(time_5, var_5, '-k', linewidth = 1, label = 'Variance')
graph_1	= ax.plot(time_6, var_6, '-k', linewidth = 1, label = 'Variance')

ax.set_xlabel('Model year')
ax.set_ylabel('Variance (Sv$^2$)')
ax.set_xlim(550, 1850)
ax.set_ylim(0, 0.8)
ax.grid()

ax2 = ax.twinx()

ax2.plot(time_4, trend_auto_corr_4 * time_4 + base_auto_corr_4, '--', color = 'dodgerblue', linewidth = 1)
ax2.plot(time_5, trend_auto_corr_5 * time_5 + base_auto_corr_5, '--', color = 'dodgerblue', linewidth = 1)
ax2.plot(time_6, trend_auto_corr_6 * time_6 + base_auto_corr_6, '--', color = 'dodgerblue', linewidth = 1)

graph_2	= ax2.plot(time_4, auto_corr_4, '-b', linewidth = 1, label = 'Auto-correlation')
graph_2	= ax2.plot(time_5, auto_corr_5, '-b', linewidth = 1, label = 'Auto-correlation')
graph_3	= ax2.plot(time_6, auto_corr_6, '-b', linewidth = 1, label = 'Auto-correlation')

ax2.set_ylabel('Lag-1 auto-correlation', color = 'b')
ax2.set_ylim(0, 1)

for tl in ax2.get_yticklabels():
    tl.set_color('b')

graphs	      = graph_1 + graph_2

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.text(0.155, 0.035, '$F_{\mathrm{ovS}} > 0$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.465, 0.035, r'$F_{\mathrm{ovS}} \approx 0$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.77, 0.035, r'$F_{\mathrm{ovS}} < 0$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)

#Significance of trends (period 1)
ax.text(0.155, 0.82, '$p < 0.01$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.155, 0.77, '$p < 0.05$', verticalalignment='center', horizontalalignment='center', color = 'b', fontsize = 12, transform = ax.transAxes)

#Significance of trends (period 2)
ax.text(0.465, 0.82, '$p < 0.05$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.465, 0.77, '$p < 0.1$', verticalalignment='center', horizontalalignment='center', color = 'b', fontsize = 12, transform = ax.transAxes)

#Significance of trends (period 3)
ax.text(0.77, 0.82, '$p < 0.05$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.77, 0.77, '$p < 0.01$', verticalalignment='center', horizontalalignment='center', color = 'b', fontsize = 12, transform = ax.transAxes)

print()
print(p_var_4, p_var_5, p_var_6)
print(p_auto_corr_4, p_auto_corr_5, p_auto_corr_6)


ax.set_title('c) AMOC strength, variance and auto-correlation')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([time_7[0], time_7[-1]], -100, 100, alpha=0.25, edgecolor='orange', facecolor='orange')
ax.fill_between([time_8[0], time_8[-1]], -100, 100, alpha=0.25, edgecolor='cyan', facecolor='cyan')
ax.fill_between([time_9[0], time_9[-1]], -100, 100, alpha=0.25, edgecolor='red', facecolor='red')

ax.plot(time_7, (trend_var_7 * time_7 + base_var_7)*10**3.0, '--', color = 'gray', linewidth = 1)
ax.plot(time_8, (trend_var_8 * time_8 + base_var_8)*10**3.0, '--', color = 'gray', linewidth = 1)
ax.plot(time_9, (trend_var_9 * time_9 + base_var_9)*10**3.0, '--', color = 'gray', linewidth = 1)

graph_1	= ax.plot(time_7, var_7*10**3.0, '-k', linewidth = 1, label = 'Variance')
graph_1	= ax.plot(time_8, var_8*10**3.0, '-k', linewidth = 1, label = 'Variance')
graph_1	= ax.plot(time_9, var_9*10**3.0, '-k', linewidth = 1, label = 'Variance')

ax.set_xlabel('Model year')
ax.set_ylabel(r'Variance ($\times 10^{-3}$ Sv$^2$)')
ax.set_xlim(550, 1850)
ax.set_ylim(0, 0.3)
ax.grid()

ax2 = ax.twinx()

ax2.plot(time_7, trend_auto_corr_7 * time_7 + base_auto_corr_7, '--', color = 'dodgerblue', linewidth = 1)
ax2.plot(time_8, trend_auto_corr_8 * time_8 + base_auto_corr_8, '--', color = 'dodgerblue', linewidth = 1)
ax2.plot(time_9, trend_auto_corr_9 * time_9 + base_auto_corr_9, '--', color = 'dodgerblue', linewidth = 1)

graph_2	= ax2.plot(time_7, auto_corr_7, '-b', linewidth = 1, label = 'Auto-correlation')
graph_2	= ax2.plot(time_8, auto_corr_8, '-b', linewidth = 1, label = 'Auto-correlation')
graph_3	= ax2.plot(time_9, auto_corr_9, '-b', linewidth = 1, label = 'Auto-correlation')

ax2.set_ylabel('Lag-1 auto-correlation', color = 'b')
ax2.set_ylim(0, 1)

for tl in ax2.get_yticklabels():
    tl.set_color('b')

graphs	      = graph_1 + graph_2

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.text(0.155, 0.035, '$F_{\mathrm{ovS}} > 0$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.465, 0.035, r'$F_{\mathrm{ovS}} \approx 0$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.77, 0.035, r'$F_{\mathrm{ovS}} < 0$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)

#Significance of trends (period 1)
ax.text(0.155, 0.81, '$p > 0.1$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.155, 0.76, '$p > 0.1$', verticalalignment='center', horizontalalignment='center', color = 'b', fontsize = 12, transform = ax.transAxes)

#Significance of trends (period 2)
ax.text(0.465, 0.81, '$p > 0.1$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.465, 0.76, '$p > 0.1$', verticalalignment='center', horizontalalignment='center', color = 'b', fontsize = 12, transform = ax.transAxes)

#Significance of trends (period 3)
ax.text(0.77, 0.95, '$p < 0.01$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 12, transform = ax.transAxes)
ax.text(0.77, 0.9, '$p > 0.1$', verticalalignment='center', horizontalalignment='center', color = 'b', fontsize = 12, transform = ax.transAxes)

print()
print(p_var_7, p_var_8, p_var_9)
print(p_auto_corr_7, p_auto_corr_8, p_auto_corr_9)

ax.set_title('d) $F_{\mathrm{ovS}}$, variance and auto-correlation')

show()

