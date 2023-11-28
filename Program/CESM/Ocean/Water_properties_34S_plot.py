#Program plots the section along 34S

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

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

year_start	= 1
year_end	= 50

#year_start	= 1701
#year_end	= 1750

#year_start	= 2151
#year_end	= 2200

depth_min 	= 0
depth_max	= 6000


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
year, lon_u, lat_u, lon_t, lat_t, depth, layer, grid_x_u, grid_x_t, v_vel_year, salt_year 	= ReadinData(files[0], depth_min_index, depth_max_index)
layer_field_u					= ma.masked_all((len(depth), len(lon_u)))
layer_field_t					= ma.masked_all((len(depth), 2, len(lon_t)))

for depth_i in range(len(depth)):
	#Determine the total length of the section, based on non-masked elements
	layer_field_u[depth_i]	= layer[depth_i]
	layer_field_u[depth_i]	= ma.masked_array(layer_field_u[depth_i], mask = v_vel_year[0, depth_i].mask)
	layer_field_t[depth_i]	= layer[depth_i]
	layer_field_t[depth_i]	= ma.masked_array(layer_field_t[depth_i], mask = salt_year[0, depth_i].mask)

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
layer_field_u_norm  = ma.masked_all(shape(layer_field_u))
grid_x_u_norm	    = ma.masked_all((len(depth), len(lon_u)))
grid_x_t_norm	    = ma.masked_all((len(depth), 2, len(lon_t)))
lat_weight	    = ma.masked_all((len(depth), 2))

for depth_i in range(len(depth)):
	#Normalise each layer
	layer_field_u_norm[depth_i]		= layer_field_u[depth_i] / np.sum(layer_field_u[depth_i])
    
	#Normalise the length
	grid_x_u_depth          		= ma.masked_array(grid_x_u, mask = v_vel_year[0, depth_i].mask)
	grid_x_u_norm[depth_i]  		= grid_x_u_depth / np.sum(grid_x_u_depth)
	grid_x_t_depth          		= ma.masked_array(grid_x_t, mask = salt_year[0, depth_i].mask)

	#Now get the lat weights
	lat_weight[depth_i]			= np.sum(grid_x_t_depth, axis = 1) / np.sum(grid_x_t_depth)

	for lat_i in range(2):
		grid_x_t_norm[depth_i, lat_i]  = grid_x_t_depth[lat_i] / np.sum(grid_x_t_depth[lat_i])

#-----------------------------------------------------------------------------------------

#Define empty array's
time_all        = np.zeros(len(files) * len(year))
transport_salt	= ma.masked_all(len(time_all))
vel_all		    = ma.masked_all((len(time_all), len(depth), 1))
vel_salt_all	= ma.masked_all((len(time_all), len(depth), 1))
salt_all	    = ma.masked_all((len(time_all), len(depth), len(lon_t)))

for file_i in range(len(files)):
    #Now determine for each year
    print(file_i)
	    
    year, lon_u, lat_u, lon_t, lat_t, depth, layer, grid_x_u, grid_x_t, v_vel_year, salt_year = ReadinData(files[file_i], depth_min_index, depth_max_index)

    for year_i in range(len(year)):
        
        #Save the time
        time_all[file_i*len(year)+year_i] = year[year_i]
        v_vel                             = v_vel_year[year_i]
        salt                              = salt_year[year_i] 

        #-----------------------------------------------------------------------------------------

    	#Determine the meridional transport
        transport	= v_vel * layer_field_u * grid_x_u

        #Determine the section averaged velocity (barotropic)
        vel_barotropic	= np.sum(transport) / np.sum(layer_field_u * grid_x_u)

    	#Determine the overturning velocity (baroclinic)
        vel_baroclinic	 = v_vel - vel_barotropic

        #Determine the zonal means
        salt_zonal      = np.sum(salt * grid_x_t_norm, axis = 2)  - 35.0
        transport_clin	= np.sum(vel_baroclinic * layer_field_u * grid_x_u, axis = 1)

        #Take the mean over the two latitudes for the salt
        salt_zonal      = np.sum(salt_zonal * lat_weight, axis = 1)

        #Save the salinity transport
        transport_salt[file_i*len(year)+year_i]	     = (-1.0 / 35.0) * np.sum(transport_clin * salt_zonal) / 1000000.0 
        vel_all[file_i*len(year)+year_i, :, 0]	     = np.sum(vel_baroclinic * grid_x_u_norm, axis = 1) * 100.0
        vel_salt_all[file_i*len(year)+year_i, :, 0] = (-1.0 / 35.0) * transport_clin * salt_zonal / 10**6.0
        salt_all[file_i*len(year)+year_i]		     = salt[:, 0]

#-----------------------------------------------------------------------------------------

time_start	    = (np.abs(time_all - year_start)).argmin()
time_end	    = (np.abs(time_all - (year_end+1))).argmin()

time_all	    = time_all[time_start:time_end]
transport_salt	= np.mean(transport_salt[time_start:time_end])
vel_all			= np.mean(vel_all[time_start:time_end, :, 0], axis = 0)
vel_salt_all	= np.mean(vel_salt_all[time_start:time_end, :, 0], axis = 0)
vel_salt_all	= vel_salt_all / layer * 1000.0
salt_all		= np.mean(salt_all[time_start:time_end], axis = 0)
#-----------------------------------------------------------------------------------------
#Get the water properties

#North Atlantic Deep Water (NADW) has negative meridional velocities
depth_index_NADW = np.where((depth >= 700) & (vel_all <= 0))[0][0]

#Antarctic bottom water (AABW) is directly below the NADW, get the first index
depth_index_AABW	= np.where((depth >= 3000) & (vel_all >= 0))[0][0]

#The Antarctic Intermediate water is between the NADW and 500 m
depth_index_AAIW	= np.where(depth >= 500)[0][0]

depth_top	= np.zeros(len(depth))

for depth_i in range(1, len(depth)):
	depth_top[depth_i]	= depth_top[depth_i - 1] + layer[depth_i - 1]

depth_AAIW	= depth_top[depth_index_AAIW]
depth_NADW	= depth_top[depth_index_NADW]
depth_AABW	= depth_top[depth_index_AABW]

lon_AAIW_index		= np.where(salt_all[depth_index_AAIW].mask == False)[0]
lon_NADW_index		= np.where(salt_all[depth_index_NADW].mask == False)[0]
lon_AABW_index		= np.where(salt_all[depth_index_AABW].mask == False)[0]
lon_AAIW_1, lon_AAIW_2	= lon_t[lon_AAIW_index[0]], lon_t[lon_AAIW_index[-1]]
lon_NADW_1, lon_NADW_2	= lon_t[lon_NADW_index[0]], lon_t[lon_NADW_index[-1]]
lon_AABW_1, lon_AABW_2	= lon_t[lon_AABW_index[0]], lon_t[lon_AABW_index[-1]]

#-----------------------------------------------------------------------------------------

depth_crop			= 1000
factor_depth_crop		= 4
depth[depth > depth_crop] 	= ((depth[depth > depth_crop] - depth_crop) / factor_depth_crop) + depth_crop

if depth_AAIW > depth_crop:
	depth_AAIW	= ((depth_AAIW - depth_crop) / factor_depth_crop) + depth_crop
if depth_NADW > depth_crop:
	depth_NADW	= ((depth_NADW - depth_crop) / factor_depth_crop) + depth_crop
if depth_AABW > depth_crop:
	depth_AABW	= ((depth_AABW - depth_crop) / factor_depth_crop) + depth_crop
#-----------------------------------------------------------------------------------------


cNorm  		= colors.Normalize(vmin=-1, vmax= 1) 		#Probablility
scalarMap 	= cm.ScalarMappable(norm=cNorm, cmap='RdBu_r') 	#Using colormap
color_south 	= scalarMap.to_rgba(-0.5)
color_north 	= scalarMap.to_rgba(0.5)

fig, ax	= subplots()

if year_start < 1750:
	ax.axhline(y = depth_AAIW, linestyle = '--', linewidth = 2.0, color = 'k')
	ax.axhline(y = depth_NADW, linestyle = '--', linewidth = 2.0, color = 'k')
	ax.axhline(y = depth_AABW, linestyle = '--', linewidth = 2.0, color = 'k')
ax.plot(vel_all, depth, '-k', linewidth = 2.0)

ax.set_xlim(-1, 1)
ax.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)
ax.grid()

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] > depth_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - depth_crop) * factor_depth_crop) + depth_crop

labels	= labels.astype(int)
ax.set_yticklabels(labels)

ax.fill_betweenx(depth, vel_all, where = vel_all >= 0.0, color = color_north, alpha = 0.50)	
ax.fill_betweenx(depth, vel_all, where = vel_all <= 0.0, color = color_south, alpha = 0.50)	

ax.set_xlabel('Meridional velocity (cm s$^{-1}$)')
ax.set_ylabel('Depth (m)')
ax.axvline(x = 0, linestyle = '--', color = 'k')

if year_start < 1750:
	ax.text(-0.95, 350, 'ASW', verticalalignment='center', horizontalalignment='left', color = 'k', fontsize=16)
	ax.text(-0.95, 750, 'AAIW', verticalalignment='center', horizontalalignment='left', color = 'k', fontsize=16)
	ax.text(-0.95, 1350, 'NADW', verticalalignment='center', horizontalalignment='left', color = 'k', fontsize=16)
	ax.text(-0.95, 1900, 'AABW', verticalalignment='center', horizontalalignment='left', color = 'k', fontsize=16)

if year_start == 1 and year_end == 50:
	ax.set_title('a) Meridional velocity (1 - 50)')

if year_start == 1701 and year_end == 1750:
	ax.set_title('b) Meridional velocity (1701 - 1750)')

if year_start == 2151 and year_end == 2200:
	ax.set_title('c) Meridional velocity (2151 - 2200)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-60, 20], y1 = np.zeros(2) + depth[0], y2 = np.zeros(2) + 2*depth[-1], color = 'gray', alpha = 0.50)

if year_start < 1750 or year_start > 4100:
	ax.plot([lon_AAIW_1, lon_AAIW_2], [depth_AAIW, depth_AAIW], linestyle = '--', linewidth = 2.0, color = 'k')
	ax.plot([lon_NADW_1, lon_NADW_2], [depth_NADW, depth_NADW], linestyle = '--', linewidth = 2.0, color = 'k')
	ax.plot([lon_AABW_1, lon_AABW_2], [depth_AABW, depth_AABW], linestyle = '--', linewidth = 2.0, color = 'k')

CS	= contourf(lon_t, depth, salt_all, levels = np.arange(34, 36.01, 0.1), extend = 'both', cmap = 'BrBG_r')
cbar	= colorbar(CS, ticks = np.arange(34, 36.01, 0.5))
cbar.set_label('Salinity (g kg$^{-1}$)')

ax.set_xlim(-60, 20)
ax.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)
ax.set_ylabel('Depth (m)')	

ax.set_xticks(np.arange(-60, 21, 10))
ax.set_xticklabels(['60$^{\circ}$W', '50$^{\circ}$W', '40$^{\circ}$W', '30$^{\circ}$W', '20$^{\circ}$W', '10$^{\circ}$W','0$^{\circ}$', '10$^{\circ}$E', '20$^{\circ}$E'])

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] > depth_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - depth_crop) * factor_depth_crop) + depth_crop

labels	= labels.astype(int)
ax.set_yticklabels(labels)

if year_start < 1750:
	ax.text(-18, 350, 'ASW', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=16)
	ax.text(-18, 750, 'AAIW', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=16)
	ax.text(-18, 1350, 'NADW', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=16)
	ax.text(-18, 1900, 'AABW', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=16)

if year_start == 1 and year_end == 50:
	ax.set_title('d) Salinity (1 - 50)')

if year_start == 1701 and year_end == 1750:
	ax.set_title('e) Salinity (1701 - 1750)')

if year_start == 2151 and year_end == 2200:
	ax.set_title('f) Salinity (2151 - 2200)')

#-----------------------------------------------------------------------------------------

cNorm  		= colors.Normalize(vmin=34, vmax= 36) 		#Probablility
scalarMap 	= cm.ScalarMappable(norm=cNorm, cmap='BrBG_r') 	#Using colormap
color_fresh 	= scalarMap.to_rgba(34.5)
color_salt 	= scalarMap.to_rgba(35.5)

fig, ax	= subplots()

if year_start < 1750:
	ax.axhline(y = depth_AAIW, linestyle = '--', linewidth = 2.0, color = 'k')
	ax.axhline(y = depth_NADW, linestyle = '--', linewidth = 2.0, color = 'k')
	ax.axhline(y = depth_AABW, linestyle = '--', linewidth = 2.0, color = 'k')
ax.plot(vel_salt_all, depth, '-k', linewidth = 2.0)

ax.set_xlim(-0.5, 0.5)
ax.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)
ax.grid()

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] > depth_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - depth_crop) * factor_depth_crop) + depth_crop

labels	= labels.astype(int)
ax.set_yticklabels(labels)

ax.set_xlabel(r'Freshwater transport (mSv m$^{-1}$)')
ax.set_ylabel('Depth (m)')
ax.axvline(x = 0, linestyle = '--', color = 'k')

ax.fill_betweenx(depth, vel_salt_all, where = vel_salt_all >= 0.0, color = color_fresh, alpha = 0.50)	
ax.fill_betweenx(depth, vel_salt_all, where = vel_salt_all <= 0.0, color = color_salt, alpha = 0.50)

if year_start < 1750:
	ax.text(-0.48, 350, 'ASW', verticalalignment='center', horizontalalignment='left', color = 'k', fontsize=16)
	ax.text(-0.48, 750, 'AAIW', verticalalignment='center', horizontalalignment='left', color = 'k', fontsize=16)
	ax.text(-0.48, 1350, 'NADW', verticalalignment='center', horizontalalignment='left', color = 'k', fontsize=16)
	ax.text(-0.48, 1900, 'AABW', verticalalignment='center', horizontalalignment='left', color = 'k', fontsize=16)	

if year_start == 1 and year_end == 50:
	ax.set_title('g) Freshwater transport (1 - 50, $F_{\mathrm{ovS}}$ = '+str(np.round(transport_salt, 2))+' Sv)')

if year_start == 1701 and year_end == 1750:
	ax.set_title('h) Freshwater transport (1701 - 1750, $F_{\mathrm{ovS}}$ = '+str(np.round(transport_salt, 2))+' Sv)')

if year_start == 2151 and year_end == 2200:
	ax.set_title('i) Freshwater transport (2151 - 2200, $F_{\mathrm{ovS}}$ = '+str(np.round(transport_salt, 2))+' Sv)')

#-----------------------------------------------------------------------------------------

show()
