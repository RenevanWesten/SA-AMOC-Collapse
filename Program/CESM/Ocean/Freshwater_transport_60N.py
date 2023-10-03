#Program determines the gyre freshwater transport

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
	#The CESM grid is structured as
	#S - S -
	#- U* - U = 34S
	#S* - S - 
	#Where the stars have the same index

	fh = netcdf.Dataset(filename, 'r')

	#First get the u-grid
	year 		= fh.variables['year'][:]							            #Model year
	lon_u 		= fh.variables['ULONG'][:]							            #Longitude
	lat_u 		= fh.variables['ULAT'][:]						                #Latitude 
	depth   	= fh.variables['z_t'][depth_min_index:depth_max_index] 			#Depth (m)
	layer		= fh.variables['dz'][depth_min_index:depth_max_index] 			#Layer thickness (m)
	grid_x_u	= fh.variables['DXU'][:] 		                               	#Zonal grid cell length (m)
	v_vel 		= fh.variables['VVEL'][:, depth_min_index:depth_max_index, 0]  #Meridional velocity (m/s)

	#Get the t-grid
	lon_t 		= fh.variables['TLONG'][:]							          #Longitude
	lat_t 		= fh.variables['TLAT'][:]					                  #Latitude 
	grid_x_t	= fh.variables['DXT'][:] 		                   	          #Zonal grid cell length (m)
	salt		= fh.variables['SALT'][:, depth_min_index:depth_max_index]   #Salinity (g / kg)
	

	fh.close()
    
	return year, lon_u, lat_u, lon_t, lat_t, depth, layer, grid_x_u, grid_x_t, v_vel, salt
    			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min 	= 0
depth_max	= 6000

#-----------------------------------------------------------------------------------------

files = glob.glob(directory+'Data/FOV_section_60N/CESM_year_*.nc')
files.sort()

#-----------------------------------------------------------------------------------------

#Get all the relevant indices to determine the mass transport
fh = netcdf.Dataset(files[0], 'r')

lat_u 		= fh.variables['ULAT'][:]	#Latitude  
depth   	= fh.variables['z_t'][:]	#Depth (m)
depth_u 	= fh.variables['HU'][:] 	#Depth at u-grid (m)
depth_t 	= fh.variables['HT'][:] 	#Depth at t-grid (m)
area		= fh.variables['TAREA'][:]	#Area at T-grid
	
fh.close()

#Get the dimensions of depth and latitude
depth_min_index 	= (fabs(depth_min - depth)).argmin()
depth_max_index 	= (fabs(depth_max - depth)).argmin() + 1

#-----------------------------------------------------------------------------------------
#Determine the section length per depth layer
year, lon_u, lat_u, lon_t, lat_t, depth, layer, grid_x_u, grid_x_t, v_vel_all, salt_all 	= ReadinData(files[0], depth_min_index, depth_max_index)
layer_field_u					= ma.masked_all((len(depth), len(lon_u[0])))
layer_field_t					= ma.masked_all((len(depth), 2, len(lon_t[0])))

for depth_i in range(len(depth)):
	#Determine the total length of the section, based on non-masked elements
	layer_field_u[depth_i]	= layer[depth_i]
	layer_field_u[depth_i]	= ma.masked_array(layer_field_u[depth_i], mask = v_vel_all[0, depth_i].mask)
	layer_field_t[depth_i]	= layer[depth_i]
	layer_field_t[depth_i]	= ma.masked_array(layer_field_t[depth_i], mask = salt_all[0, depth_i].mask)

	#Determine where the layer needs to be adjusted, partial depth cells
	depth_diff_u	= np.sum(layer_field_u, axis = 0) - depth_u
	depth_diff_t	= np.sum(layer_field_t, axis = 0) - depth_t

	#If the depth difference is negative (i.e. bottom is not reached), set to zero
	depth_diff_u	= ma.masked_where(depth_diff_u < 0, depth_diff_u)
	depth_diff_u	= depth_diff_u.filled(fill_value = 0.0)
	depth_diff_t	= ma.masked_where(depth_diff_t < 0, depth_diff_t)
	depth_diff_t	= depth_diff_t.filled(fill_value = 0.0)

	#Subtract the difference of the current layer with the difference
	layer_field_u[depth_i]	= layer_field_u[depth_i] - depth_diff_u
	layer_field_t[depth_i]	= layer_field_t[depth_i] - depth_diff_t

#Normalise layer field per layer
volume_t	        = area * layer_field_t
layer_field_u_area  = ma.masked_all((len(depth), len(lon_u[0])))
layer_field_t_area  = ma.masked_all((len(depth), 2, len(lon_t[0])))
grid_x_u_norm	    = ma.masked_all((len(depth), len(lon_u[0])))	
grid_x_t_norm	    = ma.masked_all((len(depth), 2, len(lon_t[0])))

for depth_i in range(len(depth)):
	#Normalise each layer
	layer_field_u_area[depth_i]		= layer_field_u[depth_i] * grid_x_u
	layer_field_t_area[depth_i]		= layer_field_t[depth_i] * grid_x_t
    
	#Normalise the length
	grid_x_u_depth          		= ma.masked_array(layer_field_u[depth_i] * grid_x_u, mask = v_vel_all[0, depth_i].mask)
	grid_x_t_depth          		= ma.masked_array(layer_field_t[depth_i] * grid_x_t, mask = salt_all[0, depth_i].mask)

	#Normalse per depth layer
	grid_x_u_norm[depth_i]			= grid_x_u_depth / np.sum(grid_x_u_depth)

	for lat_i in range(2):
		#Normalise for each latitude
		grid_x_t_norm[depth_i, lat_i]  = grid_x_t_depth[lat_i] / np.sum(grid_x_t_depth[lat_i])

#-----------------------------------------------------------------------------------------

#Define empty array's
time_all                = np.zeros(len(files) * len(year))
transport_salt_all		= ma.masked_all(len(time_all))

for file_i in range(len(files)):
    #Now determine for each year
    print(file_i)
    	    
    year, lon_u, lat_u, lon_t, lat_t, depth, layer, grid_x_u, grid_x_t, v_vel_all, salt_all = ReadinData(files[file_i], depth_min_index, depth_max_index)

    for year_i in range(len(year)):
        
        #Save the time
        time_all[file_i*len(year)+year_i] = year[year_i]
        v_vel                             = v_vel_all[year_i]
        salt                              = salt_all[year_i]

    	#------------------------------------------------------------------------------

        salt_prime	= ma.masked_all(np.shape(v_vel))

        for depth_i in range(len(depth)):
            for lon_i in range(len(lon_u[0])):
                #Now interpolate the salinity to the u-grid by taking the volume-averaged mean
                if v_vel[depth_i, lon_i] is ma.masked or layer_field_u[depth_i, lon_i] <= 0.0:
                    continue

                #Take the 4 T-grid cells around the current U-cell
                salt_prime_grid	= salt[depth_i, :, lon_i:lon_i+2]
                volume_t_norm	= volume_t[depth_i, :, lon_i:lon_i+2] / np.sum(volume_t[depth_i, :, lon_i:lon_i+2])

                #Take the mean over the 4 grids
                salt_prime[depth_i, lon_i]	= np.sum(salt_prime_grid*volume_t_norm) - 35.0

            #Now determine the freshwater transport
            transport_salt_all[file_i*len(year)+year_i]	= (-1.0 / 35.0) * np.sum(v_vel * salt_prime * layer_field_u_area) / 10**6.0

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_control	= ax.plot(time_all, transport_salt_all, '-k', linewidth = 1.5)

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_ylim(0, 0.6)
ax.grid()

show()

#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/Freshwater_transport_section_60N.nc', 'w')

fh.createDimension('time', len(time_all))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('FW_transport', float, ('time'), zlib=True)

fh.variables['FW_transport'].longname 	= 'Freshwater transport'

fh.variables['time'].units 		    = 'Year'
fh.variables['FW_transport'].units 		= 'Sv'

#Writing data to correct variable	
fh.variables['time'][:]     	  	    = time_all
fh.variables['FW_transport'][:] 		    = transport_salt_all

fh.close()
