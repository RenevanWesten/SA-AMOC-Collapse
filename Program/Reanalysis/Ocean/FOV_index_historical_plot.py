#Program plots the historical FOV index

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

directory_glorys	= '../../../Data/Reanalysis/GLORYS/'
directory_ora_s5	= '../../../Data/Reanalysis/ORAS5/'
directory_ora_20c	= '../../../Data/Reanalysis/ORA-20C/'
directory_ecco		= '../../../Data/Reanalysis/ECCO/'
directory_soda		= '../../../Data/Reanalysis/SODA/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	FOV		    = fh.variables['F_OV'][:]	#Fresh water

	fh.close()

	return time, FOV

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
	sig_levels 	= np.arange(50, 100, 0.1) / 100.0
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

section_name	= 'section_34S'

time_glorys, FOV_glorys		= ReadinData(directory_glorys+'Ocean/FOV_index_'+section_name+'.nc')
time_ora_s5, FOV_ora_s5	    = ReadinData(directory_ora_s5+'Ocean/FOV_index_'+section_name+'.nc')
time_ora_20c, FOV_ora_20c	= ReadinData(directory_ora_20c+'Ocean/FOV_index_'+section_name+'.nc')
time_ecco, FOV_ecco		    = ReadinData(directory_ecco+'Ocean/FOV_index_'+section_name+'.nc')
time_soda, FOV_soda		    = ReadinData(directory_soda+'Ocean/FOV_index_'+section_name+'.nc')

#-----------------------------------------------------------------------------------------

time_all	= np.arange(1900, 2023)
FOV_all		= ma.masked_all((len(time_all), 5))

for time_i in range(len(time_all)):
	#Determine the FOV mean
	if np.any(time_all[time_i] == time_glorys):
		FOV_all[time_i, 0]	= FOV_glorys[(np.abs(time_glorys - time_all[time_i])).argmin()]
	if np.any(time_all[time_i] == time_ora_s5):
		FOV_all[time_i, 1]	= FOV_ora_s5[(np.abs(time_ora_s5 - time_all[time_i])).argmin()]
	if np.any(time_all[time_i] == time_ora_20c):
		FOV_all[time_i, 2]	= FOV_ora_20c[(np.abs(time_ora_20c - time_all[time_i])).argmin()]
	if np.any(time_all[time_i] == time_ecco):
		FOV_all[time_i, 3]	= FOV_ecco[(np.abs(time_ecco - time_all[time_i])).argmin()]
	if np.any(time_all[time_i] == time_soda):
		FOV_all[time_i, 4]	= FOV_soda[(np.abs(time_soda - time_all[time_i])).argmin()]


period	= np.arange(20, 44)

trend_glorys	= ma.masked_all((len(period), 43))
trend_ora_s5	= ma.masked_all((len(period), 43))
trend_ora_20c	= ma.masked_all((len(period), 43))
trend_ecco	= ma.masked_all((len(period), 43))
trend_soda	= ma.masked_all((len(period), 43))
trend_all	= ma.masked_all((len(period), 43))

sig_glorys		= ma.masked_all((len(period), 43))
sig_ora_s5	= ma.masked_all((len(period), 43))
sig_ora_20c	= ma.masked_all((len(period), 43))
sig_ecco	= ma.masked_all((len(period), 43))
sig_soda	= ma.masked_all((len(period), 43))
sig_all		= ma.masked_all((len(period), 43))

for len_i in range(len(period)):
	#Determine the length over each period
	for year_i in range(1980, 2022-period[0]-len_i+2):
		#The starting year

		if np.any(year_i == time_glorys):
			time_index	= np.where(year_i == time_glorys)[0][0]
			time		= time_glorys[time_index:time_index+period[len_i]]
			FOV		    = FOV_glorys[time_index:time_index+period[len_i]]

			if len(time) == period[len_i]:	
				#Determine the trend over each period
				trend, error, sig		= SignificantTrend(time, FOV)
				trend_glorys[len_i, year_i-1980]	= trend*1000
				sig_glorys[len_i, year_i-1980]	= sig

		if np.any(year_i == time_ora_s5):
			time_index	= np.where(year_i == time_ora_s5)[0][0]
			time		= time_ora_s5[time_index:time_index+period[len_i]]
			FOV		= FOV_ora_s5[time_index:time_index+period[len_i]]

			if len(time) == period[len_i]:	
				#Determine the trend over each period
				trend, error, sig			= SignificantTrend(time, FOV)
				trend_ora_s5[len_i, year_i-1980]	= trend*1000
				sig_ora_s5[len_i, year_i-1980]		= sig

		if np.any(year_i == time_ora_20c):
			time_index	= np.where(year_i == time_ora_20c)[0][0]
			time		= time_ora_20c[time_index:time_index+period[len_i]]
			FOV		= FOV_ora_20c[time_index:time_index+period[len_i]]
		
			if len(time) == period[len_i]:
				#Determine the trend over each period
				trend, error, sig			= SignificantTrend(time, FOV)
				trend_ora_20c[len_i, year_i-1980]	= trend*1000
				sig_ora_20c[len_i, year_i-1980]		= sig

		if np.any(year_i == time_ecco):
			time_index	= np.where(year_i == time_ecco)[0][0]
			time		= time_ecco[time_index:time_index+period[len_i]]
			FOV		= FOV_ecco[time_index:time_index+period[len_i]]

			if len(time) == period[len_i]:	
				#Determine the trend over each period
				trend, error, sig		= SignificantTrend(time, FOV)
				trend_ecco[len_i, year_i-1980]	= trend*1000
				sig_ecco[len_i, year_i-1980]	= sig

		if np.any(year_i == time_soda):
			time_index	= np.where(year_i == time_soda)[0][0]
			time		= time_soda[time_index:time_index+period[len_i]]
			FOV		= FOV_soda[time_index:time_index+period[len_i]]

			if len(time) == period[len_i]:		
				#Determine the trend over each period
				trend, error, sig		= SignificantTrend(time, FOV)
				trend_soda[len_i, year_i-1980]	= trend*1000
				sig_soda[len_i, year_i-1980]	= sig

		if np.any(year_i == time_soda):
			time_index	= np.where(year_i == time_all)[0][0]
			time		= time_all[time_index:time_index+period[len_i]]
			FOV		= np.mean(FOV_all[time_index:time_index+period[len_i]], axis = 1)
	
			if len(time) == period[len_i]:	
				#Determine the trend over each period
				trend, error, sig		= SignificantTrend(time, FOV)
				trend_all[len_i, year_i-1980]	= trend*1000
				sig_all[len_i, year_i-1980]	= sig
#-----------------------------------------------------------------------------------------

#Now take the minimum significance as a reference
sig_glorys	= np.min(sig_glorys, axis = 1)
sig_soda	= np.min(sig_soda, axis = 1)
sig_ora_s5	= np.min(sig_ora_s5, axis = 1)
sig_ora_20c	= np.min(sig_ora_20c, axis = 1)
sig_ecco	= np.min(sig_ecco, axis = 1)
sig_all		= np.min(sig_all, axis = 1)

sig_glorys_1, sig_glorys_2, sig_glorys_3	= np.where((sig_glorys >= 0.9) & (sig_glorys < 0.95))[0], np.where((sig_glorys >= 0.95) & (sig_glorys < 0.99))[0], np.where(sig_glorys >= 0.99)[0]
sig_soda_1, sig_soda_2, sig_soda_3	= np.where((sig_soda >= 0.9) & (sig_soda < 0.95))[0], np.where((sig_soda >= 0.95) & (sig_soda < 0.99))[0], np.where(sig_soda >= 0.99)[0]
sig_ora_s5_1, sig_ora_s5_2, sig_ora_s5_3	= np.where((sig_ora_s5 >= 0.9) & (sig_ora_s5 < 0.95))[0], np.where((sig_ora_s5 >= 0.95) & (sig_ora_s5 < 0.99))[0], np.where(sig_ora_s5 >= 0.99)[0]
sig_ora_20c_1, sig_ora_20c_2, sig_ora_20c_3	= np.where((sig_ora_20c >= 0.9) & (sig_ora_20c < 0.95))[0], np.where((sig_ora_20c >= 0.95) & (sig_ora_20c < 0.99))[0], np.where(sig_ora_20c >= 0.99)[0]
sig_ecco_1, sig_ecco_2, sig_ecco_3	= np.where((sig_ecco >= 0.9) & (sig_ecco < 0.95))[0], np.where((sig_ecco >= 0.95) & (sig_ecco < 0.99))[0], np.where(sig_ecco >= 0.99)[0]
sig_all_1, sig_all_2, sig_all_3	= np.where((sig_all >= 0.9) & (sig_all < 0.95))[0], np.where((sig_all >= 0.95) & (sig_all < 0.99))[0], np.where(sig_all >= 0.99)[0]

#-----------------------------------------------------------------------------------------
fig, ax	= subplots()

graph_glorys	= plot(period, np.mean(trend_glorys, axis = 1), '-', color = 'k', linewidth = 1.5, label = 'GLORYS12V1')
graph_soda	= plot(period, np.mean(trend_soda, axis = 1), '-', color = 'r', linewidth = 1.5, label = 'SODA3.15.2')
graph_ora_s5	= plot(period, np.mean(trend_ora_s5, axis = 1), '-', color = 'b', linewidth = 1.5, label = 'ORAS5')
graph_ora_20c	= plot(period, np.mean(trend_ora_20c, axis = 1), '-', color = 'c', linewidth = 1.5, label = 'ORA-20C')
graph_ecco	= plot(period, np.mean(trend_ecco, axis = 1), '-', color = 'firebrick', linewidth =1.5, label = 'ECCO-V4r4')
graph_mean	= plot(period, np.mean(trend_all, axis = 1), '-', color = 'grey', linewidth =1.5, label = 'Mean')

ax.set_xlabel('Length of trend period (year)')
ax.set_ylabel('Freshwater transport trend (mSv yr$^{-1}$)')
ax.set_ylim(-6.5, 6.5)
ax.set_xlim(19.5, 48)
ax.set_xticks([20, 25, 30, 35, 40, 45])
ax.grid()

ax.scatter(period[sig_glorys_1], np.mean(trend_glorys, axis = 1)[sig_glorys_1], s = 30, color = 'k', marker = 'o', zorder = 10)
ax.scatter(period[sig_glorys_2], np.mean(trend_glorys, axis = 1)[sig_glorys_2], s = 30, color = 'k', marker = 's', zorder = 10)
ax.scatter(period[sig_glorys_3], np.mean(trend_glorys, axis = 1)[sig_glorys_3], s = 30, color = 'k', marker = 'D', zorder = 10)

ax.scatter(period[sig_soda_1], np.mean(trend_soda, axis = 1)[sig_soda_1], s = 30, color = 'r', marker = 'o', zorder = 10)
ax.scatter(period[sig_soda_2], np.mean(trend_soda, axis = 1)[sig_soda_2], s = 30, color = 'r', marker = 's', zorder = 10)
ax.scatter(period[sig_soda_3], np.mean(trend_soda, axis = 1)[sig_soda_3], s = 30, color = 'r', marker = 'D', zorder = 10)

ax.scatter(period[sig_ora_s5_1], np.mean(trend_ora_s5, axis = 1)[sig_ora_s5_1], s = 30, color = 'b', marker = 'o', zorder = 10)
ax.scatter(period[sig_ora_s5_2], np.mean(trend_ora_s5, axis = 1)[sig_ora_s5_2], s = 30, color = 'b', marker = 's', zorder = 10)
ax.scatter(period[sig_ora_s5_3], np.mean(trend_ora_s5, axis = 1)[sig_ora_s5_3], s = 30, color = 'b', marker = 'D', zorder = 10)

ax.scatter(period[sig_ora_20c_1], np.mean(trend_ora_20c, axis = 1)[sig_ora_20c_1], s = 30, color = 'c', marker = 'o', zorder = 10)
ax.scatter(period[sig_ora_20c_2], np.mean(trend_ora_20c, axis = 1)[sig_ora_20c_2], s = 30, color = 'c', marker = 's', zorder = 10)
ax.scatter(period[sig_ora_20c_3], np.mean(trend_ora_20c, axis = 1)[sig_ora_20c_3], s = 30, color = 'c', marker = 'D', zorder = 10)

ax.scatter(period[sig_ecco_1], np.mean(trend_ecco, axis = 1)[sig_ecco_1], s = 30, color = 'firebrick', marker = 'o', zorder = 10)
ax.scatter(period[sig_ecco_2], np.mean(trend_ecco, axis = 1)[sig_ecco_2], s = 30, color = 'firebrick', marker = 's', zorder = 10)
ax.scatter(period[sig_ecco_3], np.mean(trend_ecco, axis = 1)[sig_ecco_3], s = 30, color = 'firebrick', marker = 'D', zorder = 10)

ax.scatter(period[sig_all_1], np.mean(trend_all, axis = 1)[sig_all_1], s = 30, color = 'grey', marker = 'o', zorder = 10)
ax.scatter(period[sig_all_2], np.mean(trend_all, axis = 1)[sig_all_2], s = 30, color = 'grey', marker = 's', zorder = 10)
ax.scatter(period[sig_all_3], np.mean(trend_all, axis = 1)[sig_all_3], s = 30, color = 'grey', marker = 'D', zorder = 10)


ax.text(period[8]+0.5, np.mean(trend_glorys, axis = 1)[8]+0.1, '-5.44 mSv yr$^{-1}$ ($p < 0.1$)', verticalalignment='center', horizontalalignment='left', color = 'k', fontsize = 9)
ax.text(period[8]+2, np.mean(trend_glorys, axis = 1)[8]-0.4, '(1993 - 2020)', verticalalignment='center', horizontalalignment='left', color = 'k', fontsize = 9)

ax.text(period[-3]+2.8, np.mean(trend_soda, axis = 1)[-3]+0.7, '0.92 mSv yr$^{-1}$ ($p > 0.1$)', verticalalignment='center', horizontalalignment='center', color = 'r', fontsize = 9)
ax.text(period[-3]+2.8, np.mean(trend_soda, axis = 1)[-3]+0.2, '(1980 - 2020)', verticalalignment='center', horizontalalignment='center', color = 'r', fontsize = 9)

ax.text(period[-1], np.mean(trend_ora_s5, axis = 1)[-1]-0.75, '-1.57 mSv yr$^{-1}$ ($p < 0.01$)', verticalalignment='center', horizontalalignment='center', color = 'b', fontsize = 9)
ax.text(period[-1], np.mean(trend_ora_s5, axis = 1)[-1]-1.25, '(1980 - 2022)', verticalalignment='center', horizontalalignment='center', color = 'b', fontsize = 9)

ax.text(period[10], np.mean(trend_ora_20c, axis = 1)[10]-0.81, '-1.49 mSv yr$^{-1}$ ($p < 0.01$)', verticalalignment='center', horizontalalignment='center', color = 'c', fontsize = 9)
ax.text(period[10], np.mean(trend_ora_20c, axis = 1)[10]-1.31, '(1980 - 2009)', verticalalignment='center', horizontalalignment='center', color = 'c', fontsize = 9)

ax.text(period[6]+0.5, np.mean(trend_ecco, axis = 1)[6]+0.1, '0.13 mSv yr$^{-1}$ ($p > 0.1$)', verticalalignment='center', horizontalalignment='left', color = 'firebrick', fontsize = 9)
ax.text(period[6]+2, np.mean(trend_ecco, axis = 1)[6]-0.4, '(1992 - 2017)', verticalalignment='center', horizontalalignment='left', color = 'firebrick', fontsize = 9)


ax.text(period[-1], np.mean(trend_all, axis = 1)[-1]+0.9, '-1.20 mSv yr$^{-1}$ ($p < 0.01$)', verticalalignment='center', horizontalalignment='center', color = 'grey', fontsize = 9)
ax.text(period[-1], np.mean(trend_all, axis = 1)[-1]+0.4, '(1980 - 2022)', verticalalignment='center', horizontalalignment='center', color = 'grey', fontsize = 9)



graphs	      = graph_glorys + graph_soda + graph_ora_s5 + graph_ora_20c + graph_ecco + graph_mean

legend_labels = [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

graph_1	= ax.plot([-10], [-10], linestyle = '', color = 'k', marker = 'o', zorder = 10, label = '$p < 0.1$')
graph_2	= ax.plot([-10], [-10], linestyle = '', color = 'k', marker = 's', zorder = 10, label = '$p < 0.05$')
graph_3	= ax.plot([-10], [-10], linestyle = '', color = 'k', marker = 'D', zorder = 10, label = '$p < 0.01$')

graphs	      	= graph_1 + graph_2 + graph_3

legend_labels 	= [l.get_label() for l in graphs]
legend_2	= ax.legend(graphs, legend_labels, loc = 'upper right', ncol=1, framealpha = 1.0)
ax.add_artist(legend_1)

ax.set_title('b) Historical $F_{\mathrm{ovS}}$ trend (1980 - 2022)')

#-----------------------------------------------------------------------------------------

trend, base 	= polyfit(time_all[80:], np.mean(FOV_all[80:], axis= 1), 1)

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_glorys	= plot(time_glorys, FOV_glorys, '-', color = 'k', linewidth = 1.5, label = 'GLORYS12V1 (1/12$^{\circ}$), 1993 - 2020')
graph_soda	= plot(time_soda, FOV_soda, '-', color = 'r', linewidth = 1.5, label = 'SODA3.15.2 (1/4$^{\circ}$), 1980 - 2020')
graph_ora_s5	= plot(time_ora_s5, FOV_ora_s5, '-', color = 'b', linewidth = 1.5, label = 'ORAS5 (1/4$^{\circ}$), 1958 - 2022')
graph_ora_20c	= plot(time_ora_20c, FOV_ora_20c, '-', color = 'c', linewidth = 1.5, label = 'ORA-20C (1$^{\circ}$), 1900 - 2009')
graph_ecco	= plot(time_ecco, FOV_ecco, '-', color = 'firebrick', linewidth =1.5, label = 'ECCO-V4r4 (1$^{\circ}$), 1992 - 2017')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_ylim(-0.35, 0.35)
ax.set_xlim(1900, 2022)
ax.grid()

ax.fill_between([1900, 2025], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

graphs	      = graph_glorys + graph_soda + graph_ora_s5 + graph_ora_20c + graph_ecco

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)

#-----------------------------------------------------------------------------------------
ax2 	= fig.add_axes([0.55, 0.62, 0.33, 0.20])

ax2.fill_between([1900, 2025], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

ax2.plot(time_all[80:], time_all[80:] * trend + base, '--', color = 'gray', linewidth = 1.0)
ax2.plot(time_all, np.mean(FOV_all, axis = 1), '-', color = 'gray', linewidth = 1.5, label = '$F_{\mathrm{ovS}}$')

ax2.set_xlabel('Model year')
ax2.set_ylabel('$F_{\mathrm{ovS}}$ (Sv)')
ax2.set_ylim(-0.35, 0.02)
ax2.set_xlim(1980, 2022)
ax2.set_yticks([-0.3, -0.2, -0.1, 0])
ax2.grid()

ax2.text(0.02, 0.25, '$F_{\mathrm{ovS}}$ mean trend (1980 - 2022)', verticalalignment='center', horizontalalignment='left', color = 'k', fontsize = 9, transform = ax2.transAxes)
ax2.text(0.02, 0.071, '-1.20 mSv yr$^{-1}$ ($p < 0.01$)', verticalalignment='center', horizontalalignment='left', color = 'k', fontsize = 9, transform = ax2.transAxes)

ax2.set_title('Historical $F_{\mathrm{ovS}}$ mean')

#-----------------------------------------------------------------------------------------

ax.set_title('a) Historical $F_{\mathrm{ovS}}$')

show()

