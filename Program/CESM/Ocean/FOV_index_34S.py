#Program determines the FOV index

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
	depth   	= fh.variables['z_t'][depth_min_index:depth_max_index] 		    #Depth (m)
	layer		= fh.variables['dz'][depth_min_index:depth_max_index] 		    #Layer thickness (m)
	grid_x_u	= fh.variables['DXU'][:] 		                               #Zonal grid cell length (m)
	v_vel 		= fh.variables['VVEL'][:, depth_min_index:depth_max_index, 0]  #Meridional velocity (m/s)

	#Get the t-grid
	lon_t 		= fh.variables['TLONG'][:]							            #Longitude
	lat_t 		= fh.variables['TLAT'][:]					                    #Latitude 
	grid_x_t	= fh.variables['DXT'][:] 		                                #Zonal grid cell length (m)
	salt		= fh.variables['SALT'][:, depth_min_index:depth_max_index] 	#Salinity (g / kg)

	fh.close()
    
	return year, lon_u, lat_u, lon_t, lat_t, depth, layer, grid_x_u, grid_x_t, v_vel, salt
    			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min 	    = 0
depth_max	    = 6000

#-----------------------------------------------------------------------------------------

files = glob.glob(directory+'Data/FOV_section_34S/CESM_year_*.nc')
files.sort()

#-----------------------------------------------------------------------------------------

#Get all the relevant indices to determine the mass transport
fh = netcdf.Dataset(files[0], 'r')

lat_u 		= fh.variables['ULAT'][:]	#Latitude  
depth   	= fh.variables['z_t'][:]	#Depth (m)
depth_u 	= fh.variables['HU'][:] 	#Depth at u-grid (m)
depth_t 	= fh.variables['HT'][:] 	#Depth at t-grid (m)
	
fh.close()

#Get the dimensions of depth and latitude
depth_min_index 	= (fabs(depth_min - depth)).argmin()
depth_max_index 	= (fabs(depth_max - depth)).argmin() + 1

#-----------------------------------------------------------------------------------------
#Determine the section length per depth layer
year, lon_u, lat_u, lon_t, lat_t, depth, layer, grid_x_u, grid_x_t, v_vel_all, salt_all 	= ReadinData(files[0], depth_min_index, depth_max_index)
layer_field_u					= ma.masked_all((len(depth), len(lon_u)))
layer_field_t					= ma.masked_all((len(depth), 2, len(lon_t)))

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
layer_field_u_area  = ma.masked_all((len(depth), len(lon_u)))
layer_field_t_area  = ma.masked_all((len(depth), 2, len(lon_t)))
grid_x_t_norm	    = ma.masked_all((len(depth), 2, len(lon_t)))
lat_weight	        = ma.masked_all((len(depth), 2))

for depth_i in range(len(depth)):
	#Normalise each layer
	layer_field_u_area[depth_i]		= layer_field_u[depth_i] * grid_x_u
	layer_field_t_area[depth_i]		= layer_field_t[depth_i] * grid_x_t
    
	#Normalise the length
	grid_x_t_depth          		= ma.masked_array(grid_x_t, mask = salt_all[0, depth_i].mask)

	#Now get the lat weights
	lat_weight[depth_i]			= np.sum(grid_x_t_depth, axis = 1) / np.sum(grid_x_t_depth)

	for lat_i in range(2):
		grid_x_t_norm[depth_i, lat_i]  = grid_x_t_depth[lat_i] / np.sum(grid_x_t_depth[lat_i])

#-----------------------------------------------------------------------------------------

#Define empty array's
time_all                 = np.zeros(len(files) * len(year))
transport_all		     = ma.masked_all(len(time_all))
transport_salt_all	     = ma.masked_all(len(time_all))
transport_salt_ASW_all	 = ma.masked_all(len(time_all))
transport_salt_AAIW_all	 = ma.masked_all(len(time_all))
transport_salt_NADW_all	 = ma.masked_all(len(time_all))
transport_salt_AABW_all	 = ma.masked_all(len(time_all))

for file_i in range(len(files)):
    #Now determine for each year
    print(file_i)
	    
    year, lon_u, lat_u, lon_t, lat_t, depth, layer, grid_x_u, grid_x_t, v_vel_all, salt_all = ReadinData(files[file_i], depth_min_index, depth_max_index)

    for year_i in range(len(year)):
        
        #Save the time
        time_all[file_i*len(year)+year_i] = year[year_i]
        v_vel                             = v_vel_all[year_i]
        salt                              = salt_all[year_i]

        #Determine the meridional transport
        transport	= v_vel * layer_field_u * grid_x_u
           
        #Determine the section averaged velocity (barotropic)
        vel_barotropic	= np.sum(transport) / np.sum(layer_field_u * grid_x_u)
           
        #Determine the overturning velocity (baroclinic)
        vel_baroclinic	 = v_vel - vel_barotropic
           
        #Determine the zonal means
        salt_zonal      = np.sum(salt * grid_x_t_norm, axis = 2)  - 35.0
        transport_clin	= np.sum(vel_baroclinic * layer_field_u * grid_x_u, axis = 1)
    
    	#-----------------------------------------------------------------------------------------
    	#Get the water properties
        water_prop	= ma.masked_all((len(depth), len(lon_u)))
        water_prop_salt	= ma.masked_all((len(depth), 2, len(lon_u)))
    
    	#North Atlantic Deep Water (NADW) has negative meridional velocities
        depth_index_NADW = np.where((depth >= 700) & (transport_clin <= 0))[0][0]
    
    	#Antarctic bottom water (AABW) is directly below the NADW, get the first index
        depth_index_AABW	= np.where((depth >= 3000) & (transport_clin >= 0))[0]
    
        if len(depth_index_AABW) == 0:
    		#Always assume below 4000m depths AABW water mass (in case of no AABW detection)
            depth_index_AABW	= np.where(depth >= 4000)[0][0]
        else:
            depth_index_AABW	= depth_index_AABW[0]
    
        for depth_i in range(len(depth)):
    			
            if depth_i < depth_index_NADW:
    			#Surface water
                water_prop[depth_i]	= 1.0
    
            if depth[depth_i] >= 500 and depth_i < depth_index_NADW:
    			#Antarctic Intermediate water
                water_prop[depth_i]	= 2.0
                water_prop_salt[depth_i]= 2.0
    		
            if depth_i >= depth_index_NADW and depth_i < depth_index_AABW:
    			#North Atlantic Deep Water (NADW)
                water_prop[depth_i]	= 3.0	
                water_prop_salt[depth_i]= 3.0	
    
            if depth_i >= depth_index_AABW:
    			#The AABW is defined below the NADW 
                water_prop[depth_i]	= 4.0
                water_prop_salt[depth_i]= 4.0	
    
        water_prop	= ma.masked_array(water_prop, mask = v_vel.mask)	
    	#-----------------------------------------------------------------------------------------

    	#Determine the means over the water masses
        transport_ASW	= np.sum(ma.masked_where(water_prop != 1.0, vel_baroclinic * layer_field_u * grid_x_u), axis = 1)
        transport_AAIW	= np.sum(ma.masked_where(water_prop != 2.0, vel_baroclinic * layer_field_u * grid_x_u), axis = 1)
        transport_NADW	= np.sum(ma.masked_where(water_prop != 3.0, vel_baroclinic * layer_field_u * grid_x_u), axis = 1)
        transport_AABW	= np.sum(ma.masked_where(water_prop!= 4.0, vel_baroclinic * layer_field_u * grid_x_u), axis = 1)
    
        #Take the mean over the two latitudes for the salt
        salt_zonal      = np.sum(salt_zonal * lat_weight, axis = 1)
    
        #Determine the transport per depth layer (in Sv) and take sum to determine total transport
        transport_all[file_i*len(year)+year_i]		= np.sum(transport) / 1000000.0
    	        
    	#Determine the total salinity transport
        transport_salt_all[file_i*len(year)+year_i]	= (-1.0 / 35.0) * np.sum(transport_clin * salt_zonal) / 1000000.0 
        transport_salt_ASW_all[file_i*len(year)+year_i]	= (-1.0 / 35.0) * np.sum(transport_ASW * salt_zonal) / 1000000.0 
        transport_salt_AAIW_all[file_i*len(year)+year_i]	= (-1.0 / 35.0) * np.sum(transport_AAIW * salt_zonal) / 1000000.0 
        transport_salt_NADW_all[file_i*len(year)+year_i]	= (-1.0 / 35.0) * np.sum(transport_NADW * salt_zonal) / 1000000.0 
        transport_salt_AABW_all[file_i*len(year)+year_i]	= (-1.0 / 35.0) * np.sum(transport_AABW * salt_zonal) / 1000000.0 

#-----------------------------------------------------------------------------------------


fig, ax	= subplots()

graph_control	= ax.plot(time_all, transport_salt_all, '-k', linewidth = 1.5)

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_ylim(-0.35, 0.35)
ax.grid()

show()

#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/FOV_index_section_34S.nc', 'w')

fh.createDimension('time', len(time_all))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('Transport', float, ('time'), zlib=True)
fh.createVariable('F_OV', float, ('time'), zlib=True)
fh.createVariable('F_OV_ASW', float, ('time'), zlib=True)
fh.createVariable('F_OV_AAIW', float, ('time'), zlib=True)
fh.createVariable('F_OV_NADW', float, ('time'), zlib=True)
fh.createVariable('F_OV_AABW', float, ('time'), zlib=True)

fh.variables['Transport'].longname 	= 'Volume transport'
fh.variables['F_OV'].longname 		= 'Fresh water transport'
fh.variables['F_OV_ASW'].longname 	= 'Fresh water transport (Atlantic Surface Water)'
fh.variables['F_OV_AAIW'].longname 	= 'Fresh water transport (Antarctic Intermediate Water)'
fh.variables['F_OV_NADW'].longname 	= 'Fresh water transport (North Atlantic Deep Water)'
fh.variables['F_OV_AABW'].longname 	= 'Fresh water transport (Antarctic Bottom Water)'

fh.variables['time'].units 		    = 'Year'
fh.variables['Transport'].units 	= 'Sv'
fh.variables['F_OV'].units 		    = 'Sv'
fh.variables['F_OV_ASW'].units 		= 'Sv'
fh.variables['F_OV_AAIW'].units 	= 'Sv'
fh.variables['F_OV_NADW'].units 	= 'Sv'
fh.variables['F_OV_AABW'].units 	= 'Sv'

#Writing data to correct variable	
fh.variables['time'][:]     	  	= time_all
fh.variables['Transport'][:]    	= transport_all
fh.variables['F_OV'][:] 		    = transport_salt_all
fh.variables['F_OV_ASW'][:] 		= transport_salt_ASW_all
fh.variables['F_OV_AAIW'][:] 		= transport_salt_AAIW_all
fh.variables['F_OV_NADW'][:] 		= transport_salt_NADW_all
fh.variables['F_OV_AABW'][:] 		= transport_salt_AABW_all

fh.close()
