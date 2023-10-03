#Program plots the freshwater budget of the Atlantic

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'


#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Ocean/Freshwater_content_Atlantic.nc', 'r')

time		= fh.variables['time'][:2200]     
FW_0_1000	= fh.variables['FW_0-1000m'][:2201]    	
FW_1000_6000= fh.variables['FW_1000-6000m'][:2201]    
FW_0_6000	= fh.variables['FW_0-6000m'][:2201]    	
FW_surf		= fh.variables['FW_surface'][:2200]    				

fh.close()

#Convert the freshwater content to Sv
FW_Atlantic	= np.diff(FW_0_6000) / (86400.0 * 365.0) / 10**6.0

year_average 	= 25

time_average    	= int(len(time) / year_average)
time_average		= ma.masked_all(time_average)
FW_average     	= ma.masked_all(len(time_average))

for time_i in range(len(time_average)):
	#Loop over each time slice
	time_average[time_i]	= np.mean(time[time_i*year_average:(time_i+1)*year_average])
	FW_average[time_i] 	= np.sum(FW_Atlantic[time_i*year_average:(time_i+1)*year_average])

#-----------------------------------------------------------------------------------------

fh 			= netcdf.Dataset(directory+'Ocean/Freshwater_transport_section_34S.nc', 'r')

transport_salt_34S	= fh.variables['FW_transport'][:]		

fh.close()


fh 			= netcdf.Dataset(directory+'Ocean/Freshwater_transport_section_60N.nc', 'r')

transport_salt_60N	= fh.variables['FW_transport'][:]		

fh.close()

fh 			= netcdf.Dataset(directory+'Ocean/Freshwater_transport_Mediterranean.nc', 'r')

transport_salt_med	= fh.variables['FW_transport'][:]		

fh.close()

#Determine the freshwater transport into the Atlantic Ocean
#A positive transport means net freshwater through the section, thus 60N and the Mediterranean are subtracted
FW_transport	= transport_salt_34S - transport_salt_60N - transport_salt_med

#-----------------------------------------------------------------------------------------

#Now determine the residual term
FW_mix		= FW_Atlantic - FW_transport - FW_surf

#Determine the freshwater content w.r.t. the first 50 years
FW_0_1000	= FW_0_1000[:-1] - np.mean(FW_0_1000[0:50])
FW_1000_6000= FW_1000_6000[:-1] - np.mean(FW_1000_6000[0:50])
FW_0_6000	= FW_0_6000[:-1] - np.mean(FW_0_6000[0:50])

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

plot(time, FW_Atlantic, '-b', linewidth = 0.5)
plot(time, FW_transport, '-k', linewidth = 0.5, label = 'Transport')
plot(time, FW_surf, '-r', linewidth = 0.5, label = 'Surface')

ax.set_xlim(1, 2200)
ax.set_ylim(-1, 1)
ax.grid()

ax.set_xticks([1, 500, 1000, 1500, 2000])
ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater budget (Sv)')

ax2 = ax.twinx()

graph_4	= ax2.plot(time, FW_0_6000 / 10**14.0, '-', color = 'gray', linewidth = 2, label = '$\overline{W}$')
graph_5	= ax2.plot(time, FW_0_1000 / 10**14.0, '--', color = 'dimgray', linewidth = 2, label = r'$\overline{W}_{1000\uparrow}$')
graph_6	= ax2.plot(time, FW_1000_6000 / 10**14.0, ':', color = 'dimgray', linewidth = 2, label = r'$\overline{W}_{1000\downarrow}$')

ax2.set_ylabel(r'Freshwater content difference ($\times 10^{14}$ m$^3$)')
ax2.set_ylim(-2.5, 25)

#-----------------------------------------------------------------------------------------

graph_1	= ax.plot(time, FW_transport-100, '-k', linewidth = 2, label = r'$F_{\nabla}^b$')
graph_2	= ax.plot(time, FW_surf-100, '-r', linewidth = 2, label = '$F_{\mathrm{surf}}$')
graph_3	= ax.plot(time, FW_Atlantic-100, '-b', linewidth = 2, label = r'$\frac{\mathrm{d}\overline{W}}{\mathrm{d}t}$')

graphs	      = graph_1 + graph_2 + graph_3

legend_labels = [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc=(0.01, 0.18), ncol=1, framealpha = 1.0, numpoints = 1)

graphs	      = graph_4 + graph_5 + graph_6

legend_labels = [l.get_label() for l in graphs]
legend_2      = ax.legend(graphs, legend_labels, loc=(0.7, 0.79), ncol=1, framealpha = 1.0, numpoints = 1)

ax.add_artist(legend_1)

ax.set_title('b) Atlantic freshwater budget')


show()



