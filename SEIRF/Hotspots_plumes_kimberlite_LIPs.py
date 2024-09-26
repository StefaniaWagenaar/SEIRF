# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 12:57:09 2024

@author: 6408885
"""

"""
WRITTEN BY STEFANIA WAGENAAR

This file uses the outcome of a TRM with age error and plots the predicted 
hotspot tracks, the predicted plume motion, and the predicted reconstructed 
locations of kimberlites and LIPs. 


"""

import os.path
import pandas as pd
import numpy as np  
import pygplates
import matplotlib.pyplot as plt
from pmagpy import pmag,ipmag,pmagplotlib
from matplotlib.colors import ListedColormap, BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geoTools
import gplately

from DataLoader import *
from Functions import *

"""
This hotspot file contains selected datapoints from Doubrovine et al.
2012 and O'neill et al. 2005. Usually tracks, especially their
younger parts, have multiple dated points that are close together
in space but really close together in time, often overlapping
with age errors. This likely represents the fact that plumes can
often have active hotspots over a larger spatial range, which causes 
overestimation of plume motion. To mitigate this effect, I chose
the oldest possibility, or the possibility with the smallest 
age error, when multiple datapoints with overlapping age errors
occurred. 

"""
hs_dataframe = pd.read_csv(os.path.join('data')+'/Hotspot_tracks.csv',sep=';')
interval = 10


#%% First make the interpolated dataset of the hotspot tracks
#This is only used to plot the observed track at the intervals used in the calculations


include_chains = ['Louisville', 'Tristan', 'Reunion', 'St_Helena', 'Kerguelen', 
                  'Hawaii']

interval = 10

hotspot = hs_dataframe

age_max = max(hotspot['Age'])

age_range = np.arange(interval,130+interval,interval)

hotspot_int = []

max_ages = []
for chain in include_chains:
    age = []
    
    for i in range(1,len(hotspot)):
        if hotspot['Chain'][i] == chain:
            
            age.append(hotspot['Age'][i])
    max_ages.append(np.max(age))
            
for n,chain in enumerate(include_chains):
    for time in age_range:
        for i in range(1,len(hotspot)):
            if hotspot['Chain'][i] == chain:
                if time >= hotspot['Age'][i-1] and time < hotspot['Age'][i]:
                    diff = (time - hotspot['Age'][i-1] - hotspot['Error'][i-1])/(hotspot['Age'][i]- hotspot['Error'][i]-hotspot['Age'][i-1]- hotspot['Error'][i-1])
                    
                    min_lat = ((hotspot['Lat'][i]-hotspot['Lat'][i-1])*diff)+hotspot['Lat'][i-1]
                    min_lon = ((hotspot['Lon'][i]-hotspot['Lon'][i-1])*diff)+hotspot['Lon'][i-1]
                    
                    diff = (time - hotspot['Age'][i-1] + hotspot['Error'][i-1])/(hotspot['Age'][i] + hotspot['Error'][i]-hotspot['Age'][i-1] + hotspot['Error'][i-1])
                    
                    max_lat = ((hotspot['Lat'][i]-hotspot['Lat'][i-1])*diff)+hotspot['Lat'][i-1]
                    max_lon = ((hotspot['Lon'][i]-hotspot['Lon'][i-1])*diff)+hotspot['Lon'][i-1]
                    
                    lat = ((max_lat-min_lat)/2) + min_lat
                    lon = ((max_lon-min_lon)/2) + min_lon
                    
                    hotspot_int.append({'Chain': chain, 'Lat': lat, 'Lon': lon, 'Lat_err': ((max_lat-min_lat)/2), 'Lon_err': ((max_lon-min_lon)/2),  'Age': time, 'ref_lat': hotspot['ref_lat'][i], 'ref_lon': hotspot['ref_lon'][i], 'Plateid': hotspot['Plateid'][i]})
                
                elif time > max_ages[n] and (time-hotspot['Age'][i]) <5:
                    
                    diff = (time - hotspot['Age'][i-1] - hotspot['Error'][i-1])/(hotspot['Age'][i]- hotspot['Error'][i]-hotspot['Age'][i-1]- hotspot['Error'][i-1])
                    
                    min_lat = ((hotspot['Lat'][i]-hotspot['Lat'][i-1])*diff)+hotspot['Lat'][i-1]
                    min_lon = ((hotspot['Lon'][i]-hotspot['Lon'][i-1])*diff)+hotspot['Lon'][i-1]
                    
                    diff = (time - hotspot['Age'][i-1] + hotspot['Error'][i-1])/(hotspot['Age'][i] + hotspot['Error'][i]-hotspot['Age'][i-1] + hotspot['Error'][i-1])
                    
                    max_lat = ((hotspot['Lat'][i]-hotspot['Lat'][i-1])*diff)+hotspot['Lat'][i-1]
                    max_lon = ((hotspot['Lon'][i]-hotspot['Lon'][i-1])*diff)+hotspot['Lon'][i-1]
                    
                    lat = ((max_lat-min_lat)/2) + min_lat
                    lon = ((max_lon-min_lon)/2) + min_lon
                    
                    hotspot_int.append({'Chain': chain, 'Lat': lat, 'Lon': lon, 'Lat_err': ((max_lat-min_lat)/2), 'Lon_err': ((max_lon-min_lon)/2),  'Age': time, 'ref_lat': hotspot['ref_lat'][i], 'ref_lon': hotspot['ref_lon'][i], 'Plateid': hotspot['Plateid'][i]})
                
                    
hotspot_int = pd.DataFrame(hotspot_int)

hotspot_int.to_csv(os.path.join('data')+'/Interpolated_hotspot_tracks_'+str(interval)+'My.csv', sep =';')  


#%% Compute the predicted hotspot tracks with error
"""
The code first looks to see if the datafile already exists, and if yes, it 
checks whether each chain in include_chains is in it. If the data
has been calculated, but is wrong and you want to recalculate it, you need to 
delete it from the file first
"""

hotspots_infos = [
                  ('Reunion', pygplates.MultiPointOnSphere([(-21,56)]), 56, -21,0, 70,709,'Reunion in the ', [40,75,-25,0]),
                  ('Hawaii',pygplates.MultiPointOnSphere([(19,-155)]),-155,19,0,90,901, 'Hawaii-Emperor Bend in the ',[-140,-190,15,55]),
                  ('Louisville',pygplates.MultiPointOnSphere([(-51,-138)]),-138,-51,0,90,901,'Louisville in the ',[-130,-180,-20,-50] ),
                  ('Tristan', pygplates.MultiPointOnSphere([(-37.1,-12.25)]), -39,-11,0,140,701, 'Tristan in the ', [-30,25,-40,-10]),
                  ]

#Check if the data file already exists otherwise we create it here
file_check = os.path.exists(output_dir+'/data/Hotspots.csv') 

if file_check == False:
    print('File does not exist, making one now...')
    pole_df = pd.DataFrame(columns = ['chain','hs_lat', 'hs_lon', 'time','iteration','model',
                                      'plateid'])
    pole_df.to_csv(output_dir+'/data/Hotspots.csv')
    
if file_check:
    print('File exists, loading now...')
    pole_df = pd.read_csv(output_dir+'/data/Hotspots.csv')
    
pole_dictionary = []   
    
for j,age_error in enumerate(age_error_rot_files):
    
    
    rot_file_folder, colour, label, window, end_age, iterations,marker = age_error
    
    
    
    print(label)
    
    for info in hotspots_infos:
        
        chain, seed_point, ref_lon, ref_lat, min_age, max_age, plateid, title, extent = info
    
        #Check if the input data is missing from the file
        missing_data = False
    
        if chain in np.array(pole_df['chain']):
            pole_df = pole_df.loc[pole_df["chain"]==chain]
            print('chain')
        else:
            missing_data = True
            
        if label in list(pole_df['model']):
            pole_df = pole_df.loc[pole_df["model"]==label]
            print('model')
        else:
            missing_data = True
            
        if missing_data == False:
            print('All data was computed before')
            
        elif missing_data:
            print('There is data missing! Calculating now...')
            
            iterations_array = np.arange(0,iterations,1)
            iteration_check = [25,50,75,100,125,150,175]
            
            dictionary_points = []
            
            #We first load each rotation file and compute all rotation poles, this
            #way the rotation file has to be read only once per iteration    
            for iteration in iterations_array:
                if iteration in iteration_check:
                    print('iteration: ',iteration)
                 
                rotation_file_name = rot_file_inputdir +'/'+ rot_file_folder + '/optimised_rotation_model_' + rot_file_folder + str(iteration) + 'test.rot'
            
                rotation_features = pygplates.FeatureCollection(rotation_file_name)
            
                rotation_model = pygplates.RotationModel(rotation_file_name)
                
                # A list of times to sample the motion path - from 0 to 90Ma in 10My intervals.
                times = np.arange(min_age,max_age+10,10)
            
                # Create a motion path feature.
                motion_path_feature = pygplates.Feature.create_motion_path(
                        seed_point,
                        times,
                        valid_time=(max(times), min(times)),
                        relative_plate=plateid,
                        reconstruction_plate_id=0)
                
                # Reconstruct features to this geological time.
                reconstruction_time = min_age
                
                
                
                # Reconstruct the motion path feature to the reconstruction time.
                reconstructed_motion_paths = []
                pygplates.reconstruct(motion_path_feature, rotation_model, reconstructed_motion_paths, reconstruction_time,
                    reconstruct_type=pygplates.ReconstructType.motion_path)
                
                
                
                
                
                # Iterate over all reconstructed motion paths.
                for reconstructed_motion_path in reconstructed_motion_paths:
                
                   
                    motion_path_times = reconstructed_motion_path.get_feature().get_times()
                
                    # Iterate over the points in the motion path.
                    for point_index, point in enumerate(reconstructed_motion_path.get_motion_path()):
                
                        lat, lon = point.to_lat_lon()
                        
                
                        # The first point in the path is the oldest and the last point is the youngest.
                        # So we need to start at the last time and work our way backwards.
                        time = motion_path_times[-1-point_index]
                        
                        
                        pole_dictionary.append({'chain': chain,'hs_lat': lat, 'hs_lon': lon, 'time': time,
                                                'iteration': iteration,'model': label,'plateid': plateid})
                        
            
            
new_df = pd.DataFrame(pole_dictionary)

new_df.to_csv(output_dir+'/data/Hotspots.csv', mode='a', header=False)

#%% Plot the hotspot tracks with error

dataframe_points = pd.read_csv(output_dir+'/data/Hotspots.csv',sep=',')


hs_dataframe = pd.read_csv(os.path.join('data')+'/Interpolated_hotspot_tracks_'+str(interval)+'My.csv',sep=';')


hotspots_infos = [
                  ('Reunion', pygplates.MultiPointOnSphere([(-21,56)]), 56, -21,0, 70,709,'Reunion in the ', [40,75,-25,0]),
                  ('Hawaii',pygplates.MultiPointOnSphere([(19,-155)]),-155,19,0,90,901, 'Hawaii-Emperor Bend in the ',[-140,-190,15,55]),
                  ('Louisville',pygplates.MultiPointOnSphere([(-51,-138)]),-138,-51,0,90,901,'Louisville in the ',[-130,-180,-20,-50] ),
                  ('Tristan', pygplates.MultiPointOnSphere([(-37.1,-12.25)]), -39,-11,0,140,701, 'Tristan in the ', [-30,25,-40,-10]),
                  ]

for j,age_error in enumerate(age_error_rot_files):
    
    
    rot_file_folder, colour, label, window, end_age, iterations,marker = age_error
    
    dataframe_model = dataframe_points.loc[dataframe_points['model']==label]
    
    print(label)
    
    for info in hotspots_infos:
        
        chain, seed_point, ref_lon, ref_lat, min_age, max_age, plateid, title, extent = info
        
        dataframe_chain = dataframe_model.loc[dataframe_model['chain']==chain]
        
        print(chain)
        
        plt.rcParams['axes.facecolor']='white'
        
        lat_grid=np.arange(-90,100,10)
        lon_grid=[-180., -150., -120.,  -90.,  -60.,  -30.,    0.,   30.,   60., 90.,  120.,  150.,  180.]
        
        map_axis = plt.figure(figsize=(20,16),facecolor='white',dpi=150)
        map_axis = ipmag.make_orthographic_map(central_longitude=ref_lon, central_latitude=ref_lat,
                                               figsize=(20,16),land_color='bisque',lat_grid=lat_grid,lon_grid=lon_grid)
        gl = map_axis.gridlines(xlocs=lon_grid, ylocs=lat_grid, linewidth=1,
                          color='black', linestyle='solid',draw_labels=True)
        
        gl.xlabel_style = {'size': 40}
        gl.ylabel_style = {'size': 40}

        cmap = plt.get_cmap('viridis')
        
        # #Add in the colorbar
        # ticks = [0,20,40,60,80,100,120,130]
        # sm = plt.cm.ScalarMappable(
        #             cmap='viridis', norm=BoundaryNorm(range(0,140,10),cmap.N))
        # cb = plt.colorbar(sm, orientation = 'horizontal', ticks = ticks, shrink=0.6, pad=0.03,ax=map_axis) 
        # cb.set_label(label = 'Age (Ma)',size=24,labelpad=20)
        # cb.ax.set_xticklabels(ticks)
        # cb.ax.tick_params(labelsize=20,pad=10)
        
        etopo1.imshow(ax=map_axis, interpolation="none",alpha=0.2)
        
        ### Plot the observed trail ###
        
        
        
        lats = np.array(hs_dataframe.loc[(hs_dataframe["Chain"]==chain) & (hs_dataframe["Plateid"]==plateid),"Lat"])
        lons = np.array(hs_dataframe.loc[(hs_dataframe["Chain"]==chain) & (hs_dataframe["Plateid"]==plateid),"Lon"])
        ages = np.array(hs_dataframe.loc[(hs_dataframe["Chain"]==chain) & (hs_dataframe["Plateid"]==plateid),"Age"])
        lats_error = np.abs(np.array(hs_dataframe.loc[(hs_dataframe["Chain"]==chain) & (hs_dataframe["Plateid"]==plateid),"Lat_err"]))
        lons_error = np.abs(np.array(hs_dataframe.loc[(hs_dataframe["Chain"]==chain) & (hs_dataframe["Plateid"]==plateid),"Lon_err"]))
        
        
        
        map_axis.errorbar(lons,lats,yerr=lats_error,xerr=lons_error,fmt='none',transform=ccrs.Geodetic(), ecolor='r', capsize=2)
        
        if min_age == 0:
        
            lats = np.insert(lats,0, np.float64(hs_dataframe.loc[hs_dataframe["Chain"]==chain,"ref_lat"].iloc[0]))
            lons = np.insert(lons,0,np.float64(hs_dataframe.loc[hs_dataframe["Chain"]==chain,"ref_lon"].iloc[0]))
            ages = np.insert(ages,0,0)
        
 
        map_axis.scatter(lons,lats, transform=ccrs.Geodetic(), c='r',s=500, zorder=3, marker = 's') 
        map_axis.plot(lons, lats, transform=ccrs.Geodetic(),  color = 'r',linestyle='--', lw=10,alpha=1)
        
        
        mean_lats = [np.float64(hs_dataframe.loc[hs_dataframe["Chain"]==chain,"ref_lat"].iloc[0])]
        mean_lons = [np.float64(hs_dataframe.loc[hs_dataframe["Chain"]==chain,"ref_lon"].iloc[0])]
        
        
        times = np.arange(window,max_age,window)

        cmap2 = plt.get_cmap('viridis', 13)
        
        for n,time in enumerate(times):
            
            dataframe_time = dataframe_chain.loc[dataframe_chain['time']==time]
            
            
            lons = np.array(dataframe_time['hs_lon'])
            lats = np.array(dataframe_time['hs_lat'])
            
            
            
            #Calculate the kent distribution to find the mean pole per timestep and
            #construct a uncertainty ellipse encompassing 95% of the datapoints
            kent_distribution_95 = ipmag.kent_distribution_95(dec=lons,inc=lats) 
            
            [l1,        #Length of the major axis of the ellipse
             l2,        #Length of the minor axis of the ellipse
             bearing,   #Azimuth of the major axis relative to N
             mean_lon,  #Mean longitude
             mean_lat,  #Mean latitude
             ] = ellipse_from_kent(kent_distribution_95)


            mean_lats.append(mean_lat)
            mean_lons.append(mean_lon)
            
            
               
            #Determine zorder relative to time, this way the younger datapoints will
            #plot on top of the older datapoints
            zorder = len(times)-n
            
            #Plot fill of ellipse
            plot_pole_ellipse(map_axis,kent_distribution_95, 
                      color=cmap2(n), edgecolor=cmap2(n), marker='s', 
                      markersize=20, label='', lw=1, lower=False,filled=True, 
                      facecolor=cmap2(n), alpha=0.5, zorder=zorder+10)
            
            #Plot outline of ellipse
            plot_pole_ellipse(map_axis,kent_distribution_95, 
                      color='k', edgecolor='k', marker='s', 
                      markersize=20, label='', alpha=1, lw=2, lower=False,
                      zorder=zorder-1+10,filled=False)
            
            #Plot the mean poles per timestep
            map_axis.scatter(mean_lon,mean_lat, transform=ccrs.Geodetic(), 
                             c=cmap2(n),s=1000,edgecolor='w', marker='D', 
                             zorder=zorder+120, alpha = 0.8)
        
        #Plot the motion path as a line 
        map_axis.plot(mean_lons,mean_lats,transform=ccrs.Geodetic(),color = 'k', 
                      lw=20,label=label,marker=marker,alpha=0.8,zorder=8)
        
        map_axis.set_extent(extent)

        plt.tight_layout()



        plt.savefig(output_dir+'/figures/' +'Hotspot_track_'+chain+
                '_in_the_'+label+'.png', dpi=300) 



#%% Compute the predicted plume motion with error
"""
The code first looks to see if the datafile already exists, and if yes, it 
checks whether each chain in include_chains is in it. If the data
has been calculated, but is wrong and you want to recalculate it, you need to 
delete it from the file first
"""
hs_dataframe = pd.read_csv(os.path.join('data')+'/Hotspot_tracks.csv',sep=';')

include_chains = ['Louisville', 'Tristan', 'Reunion', 'Kerguelen' ,
                  'Hawaii' , 'New_England', 'St_Helena']


#Check if the data file already exists otherwise we create it here
file_check = os.path.exists(output_dir+'/data/Plumes.csv') 

if file_check == False:
    print('File does not exist, making one now...')
    pole_df = pd.DataFrame(columns = ['chain','lat', 'lon', 'time','model',
                                      'plateid'])
    pole_df.to_csv(output_dir+'/data/Plumes.csv')
    
if file_check:
    print('File exists, loading now...')
    pole_df = pd.read_csv(output_dir+'/data/Plumes.csv')
    
pole_dictionary = [] 
    
    
for j,age_error in enumerate(age_error_rot_files):
    
    
    rot_file_folder, colour, label, window, end_age, iterations,marker = age_error

    print(label)
    
    iterations_array = np.arange(0,iterations,1)

    dictionary = []

    chains = []
    ages = []

    it_check = [25,50,75,100,125,150,175]

    for iteration in iterations_array:
        
        if iteration in it_check:
            print(iteration)

        trail_list, output = plume_migration(hs_dataframe, label,  iteration, rot_file_inputdir, include_chains, rot_file_folder)
        trail_list = pd.DataFrame(trail_list)
        trail_list = trail_list.reset_index(drop=True)
        
        
        
        for dic in output:
            dictionary.append(dic)
        
        output = pd.DataFrame(output)
        if iteration == 0:
            
            for i in range(len(output)):
                if output['Chain'][i] not in chains:
                    chains.append(output['Chain'][i])
                    ages_array = np.array(output.loc[output['Chain']==output['Chain'][i],['Age']])
                    ages_list = []
                    for age in ages_array:
                        ages_list.append(np.float64(age))
                     
                    ages.append(ages_list)


    dataframe = pd.DataFrame(dictionary)
    dataframe = dataframe.sort_values(by=['Chain','Age'])
    
    print(dataframe)
    
    for n,chain in enumerate(chains):
        print(chain)
        
        dataframe2 = dataframe.loc[dataframe['Chain']==chain]

        for time in ages[n]:
            
            dataframe3 = dataframe2.loc[dataframe2['Age'] == time]
            
            str_age = str(np.round(time,decimals=0))
            
            
            
            if len(str_age) == 5:
            
                i = int(str_age[0:2])
                
            if len(str_age) == 4:
                i = int(str_age[0])
                
            if len(str_age) == 3:
                i = 0
                
            colour2 = cmap2(i)
    
            lons = list(dataframe3['loc_lon'])
            lats = list(dataframe3['loc_lat'])
            
            extra_lons = list(dataframe3['max_loc_lon'])
            extra_lats = list(dataframe3['max_loc_lat'])
            
            for i in range(len(extra_lons)):
                extra_lons.append(list(dataframe3['min_loc_lon'])[i])
                extra_lats.append(list(dataframe3['min_loc_lat'])[i])
            
            for i in range(len(lons)):
                
                extra_lons.append(lons[i])
                extra_lats.append(lats[i])
                
            for i in range(len(extra_lons)):
                pole_dictionary.append({'chain': chain,'lat': extra_lats[i], 'lon': extra_lons[i], 'time': time,
                                        'model': label,'plateid': plateid})
                
new_df = pd.DataFrame(pole_dictionary)

new_df.to_csv(output_dir+'/data/Plumes.csv', mode='a', header=False)         



#%% Plot the predicted plume motion with error

hotspots_infos = [
                  ('Reunion', pygplates.MultiPointOnSphere([(-21,56)]), 56, -21,0, 70,709,'Reunion in the ', [40,75,-25,0]),
                  ('Hawaii',pygplates.MultiPointOnSphere([(19,-155)]),-155,19,0,90,901, 'Hawaii-Emperor Bend in the ',[-140,-190,15,55]),
                  ('Louisville',pygplates.MultiPointOnSphere([(-51,-138)]),-138,-51,0,80,901,'Louisville in the ',[-130,-180,-20,-50] ),
                  ('Tristan', pygplates.MultiPointOnSphere([(-37.1,-12.25)]), -11,-39,0,140,701, 'Tristan in the ', [-30,25,-40,-10]),
                  ]



dataframe_points = pd.read_csv(output_dir+'/data/Plumes.csv',sep=',')

interval = 10

hs_dataframe = pd.read_csv(os.path.join('data')+'/Interpolated_hotspot_tracks_'+str(interval)+'My.csv',sep=';')

for j,age_error in enumerate(age_error_rot_files):
    
    
    rot_file_folder, colour, label, window, end_age, iterations,marker = age_error
    
    dataframe_model = dataframe_points.loc[dataframe_points['model']==label]
    
    
    print(label)
    
    for info in hotspots_infos:
        
        chain, seed_point, ref_lon, ref_lat, min_age, max_age, plateid, title, extent = info
        print(chain)
        dataframe_chain = dataframe_model.loc[dataframe_model['chain']==chain]
        
        

        
        lat_grid = np.arange(-90,95,5)
        lon_grid = np.arange(-180,185,5)
        
        map_axis = plt.figure(figsize=(20,16),facecolor='white',dpi=300)
        map_axis = ipmag.make_orthographic_map(central_longitude=ref_lon, central_latitude=ref_lat,figsize=(20,16),land_color='bisque',lat_grid=lat_grid, lon_grid = lon_grid)
        
        lat_grid=np.arange(-90,100,10)
        lon_grid=[-180., -150., -120.,  -90.,  -60.,  -30.,    0.,   30.,   60., 90.,  120.,  150.,  180.]
        
        gl = map_axis.gridlines(xlocs=lon_grid, ylocs=lat_grid, linewidth=1,
                          color='black', linestyle='solid',draw_labels=True)
        
        gl.xlabel_style = {'size': 40}
        gl.ylabel_style = {'size': 40}
        
        ### Plot the observed trail ###
        
        
        
        lats = np.array(hs_dataframe.loc[(hs_dataframe["Chain"]==chain) & (hs_dataframe["Plateid"]==plateid),"Lat"])
        lons = np.array(hs_dataframe.loc[(hs_dataframe["Chain"]==chain) & (hs_dataframe["Plateid"]==plateid),"Lon"])
        ages = np.array(hs_dataframe.loc[(hs_dataframe["Chain"]==chain) & (hs_dataframe["Plateid"]==plateid),"Age"])
        lats_error = np.abs(np.array(hs_dataframe.loc[(hs_dataframe["Chain"]==chain) & (hs_dataframe["Plateid"]==plateid),"Lat_err"]))
        lons_error = np.abs(np.array(hs_dataframe.loc[(hs_dataframe["Chain"]==chain) & (hs_dataframe["Plateid"]==plateid),"Lon_err"]))
        
        
        
        map_axis.errorbar(lons,lats,yerr=lats_error,xerr=lons_error,fmt='none',transform=ccrs.Geodetic(), ecolor='r', capsize=2)
        
        if min_age == 0:
        
            lats = np.insert(lats,0, np.float64(hs_dataframe.loc[hs_dataframe["Chain"]==chain,"ref_lat"].iloc[0]))
            lons = np.insert(lons,0,np.float64(hs_dataframe.loc[hs_dataframe["Chain"]==chain,"ref_lon"].iloc[0]))
            ages = np.insert(ages,0,0)
        
 
        map_axis.scatter(lons,lats, transform=ccrs.Geodetic(), c='r',s=500, zorder=3, marker = 's') 
        map_axis.plot(lons, lats, transform=ccrs.Geodetic(),  color = 'r',linestyle='--', lw=10,alpha=1)
        
        
        
        cmap = plt.get_cmap('viridis')

        # #Add in the colorbar
        # ticks = [0,20,40,60,80,100,120,130]
        # sm = plt.cm.ScalarMappable(
        #             cmap='viridis', norm=BoundaryNorm(range(0,140,10),cmap.N))
        # cb = plt.colorbar(sm, orientation = 'horizontal', ticks = ticks, shrink=0.6, pad=0.03,ax=map_axis) 
        # cb.set_label(label = 'Age (Ma)',size=24,labelpad=20)
        # cb.ax.set_xticklabels(ticks)
        # cb.ax.tick_params(labelsize=20,pad=10)
        
        cmap2 = plt.get_cmap('viridis', 13)
        
        etopo1.imshow(ax=map_axis, interpolation="none",alpha=0.2)
        
        times = np.arange(window,max_age,window)
        
        
        mean_lats = []
        mean_lons = []
        
        for n,time in enumerate(times):
            
            dataframe_time = dataframe_chain.loc[dataframe_chain['time']==time]

            
            
            lons = np.array(dataframe_time['lon'])
            lats = np.array(dataframe_time['lat'])
            
            
            
            #Calculate the kent distribution to find the mean pole per timestep and
            #construct a uncertainty ellipse encompassing 95% of the datapoints
            kent_distribution_95 = ipmag.kent_distribution_95(dec=lons,inc=lats) 
            
            [l1,        #Length of the major axis of the ellipse
             l2,        #Length of the minor axis of the ellipse
             bearing,   #Azimuth of the major axis relative to N
             mean_lon,  #Mean longitude
             mean_lat,  #Mean latitude
             ] = ellipse_from_kent(kent_distribution_95)


            mean_lats.append(mean_lat)
            mean_lons.append(mean_lon)
            
            
               
            #Determine zorder relative to time, this way the younger datapoints will
            #plot on top of the older datapoints
            zorder = len(times)-n
            
            # #Plot fill of ellipse
            # plot_pole_ellipse(map_axis,kent_distribution_95, 
            #           color=cmap2(n), edgecolor=cmap2(n), marker='s', 
            #           markersize=20, label='', lw=1, lower=False,filled=True, 
            #           facecolor=cmap2(n), alpha=0.3, zorder=zorder)
            
            # #Plot outline of ellipse
            # plot_pole_ellipse(map_axis,kent_distribution_95, 
            #           color=cmap2(n), edgecolor=cmap2(n), marker='s', 
            #           markersize=20, label='', alpha=1, lw=1, lower=False,
            #           zorder=zorder-1,filled=False)
            
            # #Plot the mean poles per timestep
            # map_axis.scatter(mean_lon,mean_lat, transform=ccrs.Geodetic(), 
            #                  c=cmap2(n),s=200,edgecolor=colour, marker=marker, 
            #                  zorder=zorder+120)
            
            #Plot fill of ellipse
            plot_pole_ellipse(map_axis,kent_distribution_95, 
                      color=cmap2(n), edgecolor=cmap2(n), marker='s', 
                      markersize=20, label='', lw=1, lower=False,filled=True, 
                      facecolor=cmap2(n), alpha=0.5, zorder=zorder+10)
            
            #Plot outline of ellipse
            plot_pole_ellipse(map_axis,kent_distribution_95, 
                      color='k', edgecolor='k', marker='s', 
                      markersize=20, label='', alpha=1, lw=2, lower=False,
                      zorder=zorder-1+10,filled=False)
            
            #Plot the mean poles per timestep
            map_axis.scatter(mean_lon,mean_lat, transform=ccrs.Geodetic(), 
                             c=cmap2(n),s=1000,edgecolor='w', marker='D', 
                             zorder=zorder+120, alpha = 0.8)
        
        
            
        mean_lons.insert(0,ref_lon)
        mean_lats.insert(0,ref_lat)
        
        #Plot the motion path as a line 
        map_axis.plot(mean_lons,mean_lats,transform=ccrs.Geodetic(),color = 'k', 
                      lw=20,label=label,marker=marker,alpha=0.8,zorder=8)
        
        # #Plot the motion path as a line 
        # map_axis.plot(mean_lons,mean_lats,transform=ccrs.Geodetic(),color = colour, 
        #               lw=10,label=label,marker=marker,alpha=0.7,zorder=110)
        
        
        
        map_axis.set_extent(extent)

        plt.tight_layout()



        plt.savefig(output_dir+'/figures/' +'Plume_motion_'+chain+
                '_in_the_'+label+'.png', dpi=300) 
        
        
#%% Plot one figure with the LLSVP outline and the predictions for plume migration,
#kimberlite eruption locations, and LIP eruption locations

outline_llsvp = pygplates.FeatureCollection(os.path.join('data')+'/Outline_LLSVP.gpml')

include_chains = ['Louisville', 'Tristan', 'Reunion', 'Kerguelen' ,
                  'Hawaii' , 'New_England', 'St_Helena']


dataframe_points = pd.read_csv(output_dir+'/data/Plumes.csv',sep=',')




#Define the maximum age of plotted LIPs and kimberlites
max_age = 350

#Define the plates for which to exclude the kimberlites
kimb_plates = [101]

#Set map extent
extent = [-50,90,-50,10]
#extent = None

#Choose to plot an uncertainty ellipse around the kimberlite and LIP eruption
#sites. Can make the figure a bit unreadable if there are a lot.
plot_kimb_lip_uncertainty = False

for j,age_error in enumerate(age_error_rot_files):
    
    
    
    rot_file_folder, colour, label, window, end_age, iterations,marker = age_error
    
    dataframe_kimb = pd.read_csv(output_dir+'/data/Kimberlites_rot'+label+'.csv',sep=';')
    dataframe_lips = pd.read_csv(output_dir+'/data/LIPs_rot'+label+'.csv',sep=';')
    
    dataframe_model = dataframe_points.loc[dataframe_points['model']==label]
    
    if extent:
        fig = plt.figure(figsize=(21,14),facecolor='white',dpi=300)
    else:
        fig = plt.figure(figsize=(30,20),facecolor='white',dpi=300)
        
    map_axis = plt.subplot(1,1, 1,projection=ccrs.Robinson(central_longitude = 0))
    
    map_axis.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
    map_axis.add_feature(cfeature.LAND, zorder=1, color='burlywood')
    map_axis.add_feature(cfeature.OCEAN, zorder=0, facecolor='lightcyan',alpha=0.6)
    
    lat_grid=np.arange(-90,120,30)
    #lon_grid=[-180., -150., -120.,  -90.,  -60.,  -30.,    0.,   30.,   60., 90.,  120.,  150.,  180.]
    lon_grid=[-120.,  -60.,    0.,    60.,   120.,  180.]
    gl = map_axis.gridlines(xlocs=lon_grid, ylocs=lat_grid, linewidth=1,
                      color='black', linestyle='solid',draw_labels=True)
    
    gl.xlabel_style = {'size': 30}
    gl.ylabel_style = {'size': 30}
    
    
    #We just need a random rotation file to plot the LLSVP as long as the rotation
    #of Africa is 0 at 0 Ma, it's fine
    rotation_file_name = rot_file_inputdir +'/'+ rot_file_folder + '/optimised_rotation_model_' + rot_file_folder + str(1) + 'test.rot'

    rotation_model = pygplates.RotationModel(rotation_file_name)
    
    model = gplately.PlateReconstruction(rotation_model,outline_llsvp)
    
    gplot = gplately.PlotTopologies(model,time=0)
    gplot.plot_feature(map_axis,outline_llsvp)
    
    cmap = plt.get_cmap('viridis')
    
    # #Add in the colorbar
    # ticks = list(np.arange(0,400,50))
    # sm = plt.cm.ScalarMappable(
    #             cmap='viridis', norm=BoundaryNorm(range(0,max_age+10,10),cmap.N))
    # cb = plt.colorbar(sm, orientation = 'horizontal', ticks = ticks, shrink=0.6, pad=0.07,ax=map_axis) 
    # cb.set_label(label = 'Age (Ma)',size=36,labelpad=30)
    # cb.ax.set_xticklabels(ticks)
    # cb.ax.tick_params(labelsize=30,pad=10)
    
    index = max_age/10
    
    cmap2 = plt.get_cmap('viridis', int(index))
    
    ### First plot the plume motion ###
    
    for chain in include_chains:
        
        print(chain)
        
        dataframe_chain = dataframe_model.loc[dataframe_model['chain']==chain]
        
        times = np.arange(window,max_age,window)
        
        mean_lats = []
        mean_lons = []
        
        for n,time in enumerate(times):
            
            dataframe_time = dataframe_chain.loc[dataframe_chain['time']==time]

            
            
            lons = np.array(dataframe_time['lon'])
            lats = np.array(dataframe_time['lat'])
            
            if len(lats) != 0:
            
                #Calculate the kent distribution to find the mean pole per timestep and
                #construct a uncertainty ellipse encompassing 95% of the datapoints
                kent_distribution_95 = ipmag.kent_distribution_95(dec=lons,inc=lats) 
                
                [l1,        #Length of the major axis of the ellipse
                 l2,        #Length of the minor axis of the ellipse
                 bearing,   #Azimuth of the major axis relative to N
                 mean_lon,  #Mean longitude
                 mean_lat,  #Mean latitude
                 ] = ellipse_from_kent(kent_distribution_95)
    
    
                mean_lats.append(mean_lat)
                mean_lons.append(mean_lon)
                
                
                   
                #Determine zorder relative to time, this way the younger datapoints will
                #plot on top of the older datapoints
                zorder = len(times)-n+100
                
                #Plot fill of ellipse
                plot_pole_ellipse(map_axis,kent_distribution_95, 
                          color=cmap2(n), edgecolor=cmap2(n), marker='s', 
                          markersize=20, label='', lw=1, lower=False,filled=True, 
                          facecolor=cmap2(n), alpha=0.3, zorder=zorder)
                
                #Plot outline of ellipse
                plot_pole_ellipse(map_axis,kent_distribution_95, 
                          color=cmap2(n), edgecolor=cmap2(n), marker='s', 
                          markersize=20, label='', alpha=1, lw=1, lower=False,
                          zorder=zorder-1,filled=False)
                
                #Plot the mean poles per timestep
                map_axis.scatter(mean_lon,mean_lat, transform=ccrs.Geodetic(), 
                                 c=cmap2(n),s=200,edgecolor=colour, marker=marker, 
                                 zorder=zorder+120)
        
        #Plot the motion path as a line 
        map_axis.plot(mean_lons,mean_lats,transform=ccrs.Geodetic(),color = colour, 
                      lw=10,label=label,marker=marker,alpha=0.7,zorder=110)
        
    print('Hotspots done')
        
    ### Plot the kimberlites ###
    
    
    dataframe_time = dataframe_kimb.loc[dataframe_kimb['age']<max_age]
    
    numbers = np.array(dataframe_time['number'])
    
    numbers = np.sort(numbers)
    numbers = np.unique(numbers)
    print(str(len(numbers)) + ' kimberlites to plot')

    
    for number in numbers:
      
        
        dataframe_point = dataframe_time.loc[dataframe_time['number']==number]
        
        
        
        lats = np.array(dataframe_point['lat'])
        lons = np.array(dataframe_point['lon'])
        ages = np.array(dataframe_point['age'])
        plate_id = np.array(dataframe_point['plateid'])[0]
        age = ages[0]
        
        if len(kimb_plates) > 0:
            if plate_id not in kimb_plates:
                if len(lats) != 0:

                    #Calculate the kent distribution to find the mean pole per timestep and
                    #construct a uncertainty ellipse encompassing 95% of the datapoints
                    kent_distribution_95 = ipmag.kent_distribution_95(dec=lons,inc=lats) 
                    
                    [l1,        #Length of the major axis of the ellipse
                     l2,        #Length of the minor axis of the ellipse
                     bearing,   #Azimuth of the major axis relative to N
                     mean_lon,  #Mean longitude
                     mean_lat,  #Mean latitude
                     ] = ellipse_from_kent(kent_distribution_95)
                    
                   
                    str_age = str(age)
                    
                    if len(str_age) == 5:
                    
                        n = int(str_age[0:2])
                        
                    if len(str_age) == 4:
                        n = int(str_age[0])
                        
                    colour2 = cmap2(n)
               
                    #Determine zorder relative to time, this way the younger datapoints will
                    #plot on top of the older datapoints
                    zorder = len(times)-n
                    
                    if plot_kimb_lip_uncertainty:
                    
                        #Plot fill of ellipse
                        plot_pole_ellipse(map_axis,kent_distribution_95, 
                                  color=cmap2(n), edgecolor=cmap2(n), marker='s', 
                                  markersize=20, label='', lw=1, lower=False,filled=True, 
                                  facecolor=cmap2(n), alpha=0.3, zorder=zorder)
                        
                        #Plot outline of ellipse
                        plot_pole_ellipse(map_axis,kent_distribution_95, 
                                  color=cmap2(n), edgecolor=cmap2(n), marker='s', 
                                  markersize=20, label='', alpha=1, lw=1, lower=False,
                                  zorder=zorder-1,filled=False)
                    
                    #Plot the mean poles per timestep
                    map_axis.scatter(mean_lon,mean_lat, transform=ccrs.Geodetic(), 
                                     c=cmap2(n),s=200,edgecolor=cmap2(n), marker='s', 
                                     zorder=zorder+120, alpha=0.5)
        
    print('Kimberlites done')
    
    ### Plot the LIPs ###
    
    
    
    dataframe_time = dataframe_lips.loc[dataframe_lips['age']<max_age]
    
    numbers = np.array(dataframe_time['number'])
    
    numbers = np.sort(numbers)
    numbers = np.unique(numbers)
    print(str(len(numbers)) + ' LIPs to plot')
    
    for number in numbers:
        
        
        dataframe_point = dataframe_time.loc[dataframe_time['number']==number]
        
        
        
        lats = np.array(dataframe_point['lat'])
        lons = np.array(dataframe_point['lon'])
        ages = np.array(dataframe_point['age'])
        age = ages[0]
        
        if len(lats) != 0:

            #Calculate the kent distribution to find the mean pole per timestep and
            #construct a uncertainty ellipse encompassing 95% of the datapoints
            kent_distribution_95 = ipmag.kent_distribution_95(dec=lons,inc=lats) 
            
            [l1,        #Length of the major axis of the ellipse
             l2,        #Length of the minor axis of the ellipse
             bearing,   #Azimuth of the major axis relative to N
             mean_lon,  #Mean longitude
             mean_lat,  #Mean latitude
             ] = ellipse_from_kent(kent_distribution_95)
            
           
            
            
                
            str_age = str(age)
            
            if len(str_age) == 5:
            
                n = int(str_age[0:2])
                
            if len(str_age) == 4:
                n = int(str_age[0])
                
            colour2 = cmap2(n)
       
            #Determine zorder relative to time, this way the younger datapoints will
            #plot on top of the older datapoints
            zorder = len(times)-n
            
            if plot_kimb_lip_uncertainty:
            
                #Plot fill of ellipse
                plot_pole_ellipse(map_axis,kent_distribution_95, 
                          color=cmap2(n), edgecolor=cmap2(n), marker='s', 
                          markersize=20, label='', lw=1, lower=False,filled=True, 
                          facecolor=cmap2(n), alpha=0.3, zorder=zorder)
                
                #Plot outline of ellipse
                plot_pole_ellipse(map_axis,kent_distribution_95, 
                          color=cmap2(n), edgecolor=cmap2(n), marker='s', 
                          markersize=20, label='', alpha=1, lw=1, lower=False,
                          zorder=zorder-1,filled=False)
            
            #Plot the mean poles per timestep
            map_axis.scatter(mean_lon,mean_lat, transform=ccrs.Geodetic(), 
                             c=cmap2(n),s=350,edgecolor=colour, marker='*', 
                             zorder=zorder+120)
        
    print('LIPs done')
    
    plt.tight_layout()
    
    if extent:
        map_axis.set_extent(extent)
        plt.savefig(output_dir+'/figures/' +'Kimberlites_LIPs_plumes_in_the_'+label+'_zoomed.png', dpi=300) 
    else:
        plt.savefig(output_dir+'/figures/' +'Kimberlites_LIPs_plumes_in_the_'+label+'.png', dpi=300) 
    
        