#Program determines the AMOC strength at 26N

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data
directory	= '../../../Data/CESM/'

def ReadinData(filename, depth_min_index, depth_max_index):

	fh = netcdf.Dataset(filename, 'r')

	#First get the u-grid
	year 	= fh.variables['year'][:]							            #Model year
	lon 	= fh.variables['ULONG'][:]							            #Longitude
	lat 	= fh.variables['ULAT'][:]						                #Latitude 
	depth   = fh.variables['z_t'][depth_min_index:depth_max_index] 	    #Depth (m)
	layer	= fh.variables['dz'][depth_min_index:depth_max_index] 		    #Layer thickness (m)
	grid_x	= fh.variables['DXU'][:] 		                               #Zonal grid cell length (m)
	v_vel 	= fh.variables['VVEL'][:, depth_min_index:depth_max_index, 0] #Meridional velocity (m/s)

	fh.close()
    
	return year, lon, lat, depth, layer, grid_x, v_vel
    			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min 	= 0
depth_max	= 1000

#-----------------------------------------------------------------------------------------

files = glob.glob(directory+'Data/AMOC_section_26N/CESM_year_*.nc')
files.sort()

#-----------------------------------------------------------------------------------------

#Get all the relevant indices to determine the mass transport
fh = netcdf.Dataset(files[0], 'r')

lat 		= fh.variables['ULAT'][:]	#Latitude  
depth   	= fh.variables['z_t'][:]	#Depth (m)
depth_u 	= fh.variables['HU'][:] 	#Depth at u-grid (m)
	
fh.close()

#Get the dimensions of depth and latitude
depth_min_index 	= (fabs(depth_min - depth)).argmin()
depth_max_index 	= (fabs(depth_max - depth)).argmin() + 1

#-----------------------------------------------------------------------------------------
#Determine the section length per depth layer
year, lon, lat, depth, layer, grid_x, v_vel_all	= ReadinData(files[0], depth_min_index, depth_max_index)
layer_field				= ma.masked_all((len(depth), len(lon[0])))
depth_top				= np.zeros(len(depth))

for depth_i in range(1, len(depth)):
	#Determine the depth top
	depth_top[depth_i]	= depth_top[depth_i-1] + layer[depth_i-1]

for depth_i in range(len(depth)):
	#Determine the total length of the section, based on non-masked elements
	layer_field[depth_i]	= layer[depth_i]
	layer_field[depth_i]	= ma.masked_array(layer_field[depth_i], mask = v_vel_all[0, depth_i].mask)

	#Determine where the layer needs to be adjusted, partial depth cells
	depth_diff	= np.sum(layer_field, axis = 0) - depth_u

	if depth_i == len(depth) - 1:
		#Last layer, get the depth difference with respect to top and depth max boundary
		depth_diff	= layer_field[depth_i] - (depth_max -  depth_top[depth_i])

	#If the depth difference is negative (i.e. bottom is not reached), set to zero
	depth_diff	= ma.masked_where(depth_diff < 0, depth_diff)
	depth_diff	= depth_diff.filled(fill_value = 0.0)

	#Subtract the difference of the current layer with the difference
	layer_field[depth_i]	= layer_field[depth_i] - depth_diff

#-----------------------------------------------------------------------------------------
time_all            = np.zeros(len(files) * len(year))
transport_all		= ma.masked_all(len(time_all))

for file_i in range(len(files)):
    #Now determine for AMOC strength
    print(file_i)
	    
    year, lon, lat, depth, layer, grid_x, v_vel_all = ReadinData(files[file_i], depth_min_index, depth_max_index)

    for year_i in range(len(year)):
        #Determine the meridional transport
        time_all[file_i*len(year)+year_i] = year[year_i]
        transport	                      = v_vel_all[year_i] * layer_field * grid_x

        #Determine the transport per depth layer (in Sv) and take sum to determine total transport
        transport_all[file_i*len(year)+year_i]	= np.sum(transport) / 1000000.0
		        
#-----------------------------------------------------------------------------------------


fig, ax	= subplots()

graph_control	= ax.plot(time_all, transport_all, '-k', linewidth = 1.5)

ax.set_xlabel('Model year')
ax.set_ylabel('Volume transport (Sv)')
ax.set_ylim(-2, 22)
ax.grid()

show()

#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc', 'w')

fh.createDimension('time', len(time_all))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('Transport', float, ('time'), zlib=True)

fh.variables['Transport'].longname 	= 'Volume transport'

fh.variables['time'].units 		    = 'Year'
fh.variables['Transport'].units 	= 'Sv'

#Writing data to correct variable	
fh.variables['time'][:]     	  	= time_all
fh.variables['Transport'][:]    	= transport_all

fh.close()
