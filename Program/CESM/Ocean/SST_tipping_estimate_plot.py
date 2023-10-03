#Program estimates the tipping point from the monthly-averaged SST of the subpolar gyre

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats
from scipy.optimize import minimize

#Making pathway to folder with all data
directory	    = '../../../Data/CESM/'
directory_obs	= '../../../Data/HadISST/'


def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	temp_subpolar	= fh.variables['temp_subpolar'][:]	#AMOC strength (Sv)

	fh.close()

	return time, temp_subpolar

def estimate_OU(data, delta):
	"""Estimate the parameters"""

	data_up		= data[1:]
	data_low	= data[:-1]
	data_mean_0	= np.mean(data)
	alpha_0		= -np.log(np.sum((data_up - data_mean_0) * (data_low - data_mean_0)) / np.sum((data_low - data_mean_0)**2.0)) / delta
	sigma_20	= np.mean(np.diff(temp_subpolar_0)**2.0) / delta

	#Estimate the parameters from fit
	fun_loglik_OU			= lambda pars: LogLik_OU(pars, temp_subpolar_0, delta)
	res				= minimize(fun_loglik_OU, np.array([alpha_0, data_mean_0, sigma_20]))
	alpha_0, data_mean_0, sigma_20	= res.x[0], res.x[1], res.x[2]

	return alpha_0, data_mean_0, sigma_20

def LogLik_OU(pars, data, delta):
	"""Estimated parameters are alpha, mu and sigma"""

	alpha_0		= np.max([pars[0], 0.001])	#Ensuring alpha is positive
	mu_0		= pars[1]			#Mean
	sigma_2		= np.max([0, pars[2]])		#Infinitisimal variance, should be positive
	data_up		= data[1:]
	data_low	= data[:-1]
	time_2		= np.arange(1, len(data)) * delta
	gamma_2		= sigma_2 / (2 * alpha_0)
	rho_0		= np.exp(-alpha_0 * delta)	#Aut-correlation
	m_part		= data_up - data_low * rho_0 - mu_0 * (1 - rho_0)	#Observation minus the mean of transition
	v_part 		= gamma_2 * (1 - rho_0**2.0)	#Variance of transition distribution
	loglik		= - len(data) * (np.log(v_part)) - np.sum(m_part**2.0 / v_part)

	return -loglik

def LogLik(pars, data, delta, alpha_0, data_mean_0, sigma_20, pen = 0, final_fit = False):
	"""Estimated parameters are tau and a"""

	tau		= pars[0]
	a		= np.max([0.1, pars[1]])		#Drift
	m		= data_mean_0 - alpha_0 / (2 * a) 	#Constant mean shift
	lambda_0	= -alpha_0**2.0 / (4 * a)		#Stationary level of control parameter
	sigma_2		= np.copy(sigma_20)			#Infinitisimal variance
	data_up		= data[1:]
	data_low	= data[:-1]
	time_2		= np.arange(1, len(temp_subpolar_2)) * delta
	lam_seq		= lambda_0 * (1.0 - (time_2 / tau))
	alpha_seq	= 2 * np.sqrt(-a * lam_seq)
	gamma_2_seq	= sigma_2 / (2 * alpha_seq)
	rho_seq		= np.exp(-alpha_seq * delta)
	mu_seq		= m + np.sqrt(-lam_seq/a)

	if np.all(alpha_seq == alpha_seq) == False:
		return 10000

	if final_fit == True:
		return alpha_seq, gamma_2_seq, rho_seq, mu_seq

	#Calculating the Strang splitting scheme pseudo likelihood
	fh_half_tmp	= a * delta * (data_low - mu_seq) / 2.0
	fh_half		= (mu_seq * fh_half_tmp + data_low) / (fh_half_tmp + 1)
	fh_half_inv	= (mu_seq * fh_half_tmp - data_up) / (fh_half_tmp - 1)
	mu_h		= fh_half * rho_seq + mu_seq * (1.0 - rho_seq)
	m_part		= fh_half_inv - mu_h
	var_part	= gamma_2_seq * (1.0 - rho_seq**2.0)
	det_Dfh_half_inv= 1.0 / (a * delta * (data_up - mu_seq)/2.0 - 1)**2.0
	loglik 		= -np.sum(np.log(var_part)) - np.sum(m_part**2.0 / var_part) + 2*np.sum(np.log(det_Dfh_half_inv))

	if a < 1:
		loglik	= loglik - pen * len(data) * (1/a - 1)

	return -loglik


#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

window		    = 50 	#Sliding window (years)

#-----------------------------------------------------------------------------------------


data		    = np.loadtxt(directory_obs+'Ocean/AMOC_data.txt', delimiter=' ')

time 		    = data[:, 0]
temp_subpolar	= data[:, 3] #Subtract the global mean twice (as in Ditlevsen)


tau_stationary	= 1924
delta		    = 1/12.

#Stationary part
temp_subpolar_0			= temp_subpolar[np.where(time <= tau_stationary)[0]]
alpha_0, data_mean_0, sigma_20	= estimate_OU(temp_subpolar_0, delta)

#Now determine the ramping
temp_subpolar_2	        = temp_subpolar[np.where(time > tau_stationary)[0]]
fun_loglik	            = lambda pars: LogLik(pars, temp_subpolar_2, delta, alpha_0, data_mean_0, sigma_20, pen = 0.004)
res		                = minimize(fun_loglik, np.array([100, 1]), method = 'nelder-mead')

tau	        = res.x[0]
a	        = res.x[1]
m	        = data_mean_0 - alpha_0 / (2 * a)
lambda_0    = -alpha_0**2.0 / (4 * a)
tc	        = tau + tau_stationary

time_2		= np.arange(1, tau * 1/delta) * delta
lam_seq		= lambda_0 * (1.0 - (time_2 / tau))
alpha_seq	= 2 * np.sqrt(-a * lam_seq)
gamma_2_seq	= sigma_20 / (2 * alpha_seq)
rho_seq		= np.exp(-alpha_seq * delta)
mu_seq		= m + np.sqrt(-lam_seq/a)
sigma_seq	= 2 * alpha_seq * gamma_2_seq

print('Check results from Ditlevsen (2023)')
print(alpha_0, data_mean_0, sigma_20, tau, a, m, lambda_0, tc)
print()
print('Projected AMOC tipping from observed SST subpolar (monthly):', tc)
print()

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

year_start	    = 1427
year_end	    = year_start + 150
tau_stationary	= year_start + 53

fh              = netcdf.Dataset(directory+'Ocean/SST_subpolar.nc', 'r')

time		    = fh.variables['time'][(year_start-1)*12:(year_end)*12]     	  	
temp_subpolar	= fh.variables['SST_sub'][(year_start-1)*12:(year_end)*12] 		
temp_global	    = fh.variables['SST_global'][(year_start-1)*12:(year_end)*12]		

fh.close()

for month_i in range(12):
	#Remove the monthly mean
	month_index			        = np.arange(month_i, len(time), 12)
	temp_subpolar_mean		    = np.mean(temp_subpolar[month_index])
	temp_global_mean		    = np.mean(temp_global[month_index])

	temp_subpolar[month_index]	= temp_subpolar[month_index] - temp_subpolar_mean
	temp_global[month_index]	= temp_global[month_index] - temp_global_mean

temp_subpolar	= temp_subpolar - temp_global
delta		    = 1/12.

#Stationary part
temp_subpolar_0			        = temp_subpolar[np.where(time <= tau_stationary)[0]]
alpha_0, data_mean_0, sigma_20	= estimate_OU(temp_subpolar_0, delta)

#Now determine the ramping
temp_subpolar_2	= temp_subpolar[np.where(time > tau_stationary)[0]]
fun_loglik	    = lambda pars: LogLik(pars, temp_subpolar_2, delta, alpha_0, data_mean_0, sigma_20, pen = 0.004)
res		        = minimize(fun_loglik, np.array([100, 1]), method = 'nelder-mead')

tau	= res.x[0]
a	= res.x[1]
m	= data_mean_0 - alpha_0 / (2 * a)
lambda_0= -alpha_0**2.0 / (4 * a)
tc	= tau + tau_stationary

print('Projected AMOC tipping from the CESM (starting year 1427):', tc)
print()

#Now generate the fits
time_fit_1	= np.arange(1, tau * 1/delta) * delta
lam_seq		= lambda_0 * (1.0 - (time_fit_1 / tau))
alpha_seq	= 2 * np.sqrt(-a * lam_seq)
gamma_2_fit_1	= sigma_20 / (2 * alpha_seq)
rho_fit_1	= np.exp(-alpha_seq * delta)
mean_fit_1	= m + np.sqrt(-lam_seq/a)
time_fit_1	= time_fit_1 + tau_stationary

#Add the stationary part
time_fit_1	= np.insert(time_fit_1, 0, time[0])
mean_fit_1	= np.insert(mean_fit_1, 0, mean_fit_1[0])
gamma_2_fit_1	= np.insert(gamma_2_fit_1, 0, gamma_2_fit_1[0])
rho_fit_1	= np.insert(rho_fit_1, 0, rho_fit_1[0])
#-----------------------------------------------------------------------------------------

time_window_1		= np.zeros(len(time) - window * 12 + 1)
var_window_1		= np.zeros(len(time_window_1))
auto_corr_window_1	= np.zeros(len(time_window_1))

for time_i in range(len(time) - window * 12 + 1):
	#Determine the variance and auto-correlation over sliding window
	time_window_1[time_i]	= time[time_i:time_i+window*12][window * 6]

	#Remove trend to make data stationary
	data_window		= temp_subpolar[time_i:time_i+window*12]
	trend, base		= polyfit(np.arange(window*12), data_window, 1)
	data_window		= data_window - ((trend * np.arange(window*12)) + base)

	#Determine the variance
	var_window_1[time_i]		= np.var(data_window)
	auto_corr_window_1[time_i]	=  np.corrcoef(data_window[:-1], data_window[1:])[0,1]

#-----------------------------------------------------------------------------------------


fig, ax	= subplots()
plot(time, temp_subpolar, '-r', linewidth = 0.25)
plot(time_fit_1, mean_fit_1, '--', linewidth = 2.0, color = 'firebrick')

#---------------------------------------------------------------

year_start	    = 1503
year_end	    = year_start + 150
tau_stationary	= year_start + 64

fh = netcdf.Dataset(directory+'Ocean/SST_subpolar.nc', 'r')

time		= fh.variables['time'][(year_start-1)*12:(year_end)*12]     	  	
temp_subpolar	= fh.variables['SST_sub'][(year_start-1)*12:(year_end)*12] 		
temp_global	= fh.variables['SST_global'][(year_start-1)*12:(year_end)*12]		

fh.close()

for month_i in range(12):
	#Remove the monthly mean
	month_index			= np.arange(month_i, len(time), 12)
	temp_subpolar_mean		= np.mean(temp_subpolar[month_index])
	temp_global_mean		= np.mean(temp_global[month_index])

	temp_subpolar[month_index]	= temp_subpolar[month_index] - temp_subpolar_mean
	temp_global[month_index]	= temp_global[month_index] - temp_global_mean

temp_subpolar	= temp_subpolar - temp_global
delta		= 1/12.

#Stationary part
temp_subpolar_0			= temp_subpolar[np.where(time <= tau_stationary)[0]]
alpha_0, data_mean_0, sigma_20	= estimate_OU(temp_subpolar_0, delta)

#Now determine the ramping
temp_subpolar_2	= temp_subpolar[np.where(time > tau_stationary)[0]]
fun_loglik	= lambda pars: LogLik(pars, temp_subpolar_2, delta, alpha_0, data_mean_0, sigma_20, pen = 0.004)
res		= minimize(fun_loglik, np.array([100, 1]), method = 'nelder-mead')

tau	        = res.x[0]
a	        = res.x[1]
m	        = data_mean_0 - alpha_0 / (2 * a)
lambda_0    = -alpha_0**2.0 / (4 * a)
tc	        = tau + tau_stationary

print('Projected AMOC tipping from the CESM (starting year 1503):', tc)
print()

#Now generate the fits
time_fit_2	= np.arange(1, tau * 1/delta) * delta
lam_seq		= lambda_0 * (1.0 - (time_fit_2 / tau))
alpha_seq	= 2 * np.sqrt(-a * lam_seq)
gamma_2_fit_2	= sigma_20 / (2 * alpha_seq)
rho_fit_2	= np.exp(-alpha_seq * delta)
mean_fit_2	= m + np.sqrt(-lam_seq/a)
time_fit_2	= time_fit_2 + tau_stationary

#Add the stationary part
time_fit_2	= np.insert(time_fit_2, 0, time[0])
mean_fit_2	= np.insert(mean_fit_2, 0, mean_fit_2[0])
gamma_2_fit_2	= np.insert(gamma_2_fit_2, 0, gamma_2_fit_2[0])
rho_fit_2	= np.insert(rho_fit_2, 0, rho_fit_2[0])
#-----------------------------------------------------------------------------------------

time_window_2		= np.zeros(len(time) - window * 12 + 1)
var_window_2		= np.zeros(len(time_window_2))
auto_corr_window_2	= np.zeros(len(time_window_2))

for time_i in range(len(time) - window * 12 + 1):
	#Determine the variance and auto-correlation over sliding window
	time_window_2[time_i]	= time[time_i:time_i+window*12][window * 6]

	#Remove trend to make data stationary
	data_window		= temp_subpolar[time_i:time_i+window*12]
	trend, base		= polyfit(np.arange(window*12), data_window, 1)
	data_window		= data_window - ((trend * np.arange(window*12)) + base)

	#Determine the variance
	var_window_2[time_i]		= np.var(data_window)
	auto_corr_window_2[time_i]	= np.corrcoef(data_window[:-1], data_window[1:])[0,1]


#---------------------------------------------------------------

graph_1		= plot([-100, -100], [-100, 100], '-r', linewidth = 2.0, label = '1427 - 1577')
graph_2		= plot([-100, -100], [-100, 100], '-b', linewidth = 2.0, label = '1503 - 1653')
graph_3		= plot([-100, -100], [-100, 100], '--k', linewidth = 2.0, label = 'Fit')

graphs	     	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]

plot(time, temp_subpolar, '-b', linewidth = 0.25)
plot(time_fit_2, mean_fit_2, '--', linewidth = 2.0, color = 'darkblue')

ax.set_xlim(1400, 1750)
ax.set_ylim(-1.5, 1.5)
ax.set_xlabel('Model year')
ax.set_ylabel('Temperature anomaly ($^{\circ}$C)')
ax.grid()

ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('a) SST subpolar region (monthly)')
#---------------------------------------------------------------
fig, ax	= subplots()

plot(time_window_1, var_window_1, '-r')
plot(time_fit_1, gamma_2_fit_1, '--', linewidth = 2.0, color = 'firebrick')
plot(time_window_2, var_window_2, '-b')
plot(time_fit_2, gamma_2_fit_2, '--', linewidth = 2.0, color = 'darkblue')

ax.set_xlim(1400, 1750)
ax.set_ylim(0, 0.3)
ax.set_xlabel('Model year')
ax.set_ylabel('Variance ($^{\circ}$C$^2$)')
ax.grid()

ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('b) SST subpolar region (monthly), variance')

#---------------------------------------------------------------

fig, ax	= subplots()

plot(time_window_1, auto_corr_window_1, '-r')
plot(time_fit_1, rho_fit_1, '--', linewidth = 2.0, color = 'firebrick')
plot(time_window_2, auto_corr_window_2, '-b')
plot(time_fit_2, rho_fit_2, '--', linewidth = 2.0, color = 'darkblue')

ax.set_xlim(1400, 1750)
ax.set_ylim(0.8, 1)
ax.set_xlabel('Model year')
ax.set_ylabel('Lag-1 auto-correlation')
ax.set_yticks([0.8, 0.85, 0.9, 0.95, 1.0])
ax.grid()

ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('c) SST subpolar region (monthly), auto-correlation')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Now determing the estimate for all possible starting years
tipping_estimate	= ma.masked_all((151, 31))

for year_start in range(1400, 1551):
	for t_0_i in range(50, 81):

		year_end	= year_start + 150
		tau_stationary	= year_start + t_0_i

		#Read in the data
		fh = netcdf.Dataset(directory+'Ocean/SST_subpolar.nc', 'r')

		time		    = fh.variables['time'][(year_start-1)*12:(year_end)*12]     	  	
		temp_subpolar	= fh.variables['SST_sub'][(year_start-1)*12:(year_end)*12] 		
		temp_global	    = fh.variables['SST_global'][(year_start-1)*12:(year_end)*12]		

		fh.close()

		for month_i in range(12):
			#Remove the monthly mean
			month_index			= np.arange(month_i, len(time), 12)
			temp_subpolar_mean		= np.mean(temp_subpolar[month_index])
			temp_global_mean		= np.mean(temp_global[month_index])

			temp_subpolar[month_index]	= temp_subpolar[month_index] - temp_subpolar_mean
			temp_global[month_index]	= temp_global[month_index] - temp_global_mean

		temp_subpolar	= temp_subpolar - temp_global
		delta		    = 1/12.

		#Stationary part
		temp_subpolar_0			= temp_subpolar[np.where(time <= tau_stationary)[0]]
		alpha_0, data_mean_0, sigma_20	= estimate_OU(temp_subpolar_0, delta)

		#Now determine the ramping
		temp_subpolar_2	= temp_subpolar[np.where(time > tau_stationary)[0]]
		fun_loglik	= lambda pars: LogLik(pars, temp_subpolar_2, delta, alpha_0, data_mean_0, sigma_20, pen = 0.004)
		res		= minimize(fun_loglik, np.array([100, 1]), method = 'nelder-mead')

		tau	        = res.x[0]
		a	        = res.x[1]
		m	        = data_mean_0 - alpha_0 / (2 * a)
		lambda_0    = -alpha_0**2.0 / (4 * a)
		tc	        = tau + tau_stationary

		time_2		= np.arange(1, tau * 1/delta) * delta
		lam_seq		= lambda_0 * (1.0 - (time_2 / tau))
		alpha_seq	= 2 * np.sqrt(-a * lam_seq)
		gamma_2_seq	= sigma_20 / (2 * alpha_seq)
		rho_seq		= np.exp(-alpha_seq * delta)
		mu_seq		= m + np.sqrt(-lam_seq/a)
		sigma_seq	= 2 * alpha_seq * gamma_2_seq

		#Save the tipping point
		tipping_estimate[year_start - 1400, t_0_i-50]	= tc

	print('Starting year', year_start, 'estimate:', np.mean(tipping_estimate[year_start - 1400], axis = 0))

#-----------------------------------------------------------------------------------------

time	= np.arange(1400, 1551)

fig, ax	= subplots()

ax.fill_between(time, y1 = np.percentile(tipping_estimate, 2.5, axis = 1), y2 = np.percentile(tipping_estimate, 97.5, axis = 1), alpha=0.25, edgecolor='k', facecolor='k')
graph_1 = ax.plot(time, np.mean(tipping_estimate, axis = 1), '-k', linewidth = 2.0, label = 'Tipping point estimate, $t_c$')
graph_2	= ax.plot([1400, 2100], [1732.38522, 1732.35822], ':', color = 'firebrick', linewidth = 2.0, label = '$\mathrm{d}_t F_{\mathrm{ovS}} = 0$', zorder = 10)

ax.set_xlim(1400, 1550)
ax.set_ylim(1600, 2000)
ax.set_xlabel('Model year start')
ax.set_ylabel('Model year estimate to tipping point')
ax.set_xticks([1400, 1425, 1450, 1475, 1500, 1525, 1550])
ax.grid()

graphs	      = graph_1 + graph_2

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.quiver(1427, 1630, 0, 3, scale = 30, color = 'r')
ax.text(1427, 1625, '1427 - 1557', verticalalignment='top', horizontalalignment='center', color = 'r', fontsize = 12)

ax.quiver(1503, 1630, 0, 3, scale = 30, color = 'b')
ax.text(1503, 1625, '1503 - 1653', verticalalignment='top', horizontalalignment='center', color = 'b', fontsize = 12)

ax.set_title('d) Tipping point estimate from SST subpolar region')

show()
