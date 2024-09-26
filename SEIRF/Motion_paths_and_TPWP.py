# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:34:30 2024

@author: 6408885
"""

"""
WRITTEN BY STEFANIA WAGENAAR

This file uses the outcome of a TRM with age error and plots motion paths and
true polar wander paths with the associated uncertainty. When using all the data
provided in the zip file, it should be possible to not run the data loading files
and just load all input and run the cells that make figures. 


"""



import pygplates
import matplotlib.pyplot as plt
import numpy as np
import geoTools
import cartopy.crs as ccrs
import cartopy.feature
import pandas as pd
from matplotlib.colors import ListedColormap, BoundaryNorm
from pmagpy import pmag,ipmag,pmagplotlib
import rot_functions as rf
import optimised_rotation_updater

from DataLoader import *
from Functions import *


##############################################################################
## Plot the motion path of a point with the age error
##############################################################################
#%% Motion path with age error initial input and functions
 



plt.rcParams['axes.facecolor']='white'

#Coordinates of the point, default is set in Central Africa
lat_pole0 = list([-2.31])
lon_pole0 = list([19.15])

#Coordinates of the point, this option is set on the south pole
# lat_pole0 = list([-90])
# lon_pole0 = list([0])

#See if latitude is -90, because if so, the angle of the rotation pole has to 
#be inverted
test = False
if lat_pole0[0] == -90:
    test = True

#Plate id of the plate where the point is located
plate_id = 701





#%% Motion path with age error load data
#### First load the data ####
"""
The code first looks to see if the datafile already exists, and if yes, it 
checks whether the data put in the previous cell is already in it. If the data
has been calculated, but is wrong and you want to recalculate it, you need to 
delete it from the file first
"""

#Check if the data file already exists otherwise we create it here
file_check = os.path.exists(output_dir+'/data/Motion_paths.csv') 

if file_check == False:
    print('File does not exist, making one now...')
    pole_df = pd.DataFrame(columns = ['lat', 'lon', 'time','iteration','model',
                                      'plateid','actual_time','lat_pole','lon_pole'])
    pole_df.to_csv(output_dir+'/data/Motion_paths.csv')
    
if file_check:
    print('File exists, loading now...')
    pole_df = pd.read_csv(output_dir+'/data/Motion_paths.csv')
    

    

    

pole_dictionary = []

for j,age_error in enumerate(age_error_rot_files):
    
    
    rot_file_folder, colour, label, window, end_age, iterations,marker = age_error
    
    
    
    print(label)
    
    #Check if the input data is missing from the file
    missing_data = False
    
    
        
    if plate_id in np.array(pole_df['plateid']):
        pole_df = pole_df.loc[pole_df["plateid"]==plate_id]
        print('plateid')
    else:
        missing_data = True
        
    if label in list(pole_df['model']):
        pole_df = pole_df.loc[pole_df["model"]==label]
        print('model')
    else:
        missing_data = True
        
    if lat_pole0[0] in np.array(pole_df['lat_pole']):
        pole_df = pole_df.loc[pole_df["lat_pole"]==lat_pole0[0]]
        print('lat_pole0')
    else:
        missing_data = True
    
    if missing_data == False:
        print('All data was computed before')
        
    elif missing_data:
        print('There is data missing! Calculating now...')
        
        #Can be used to plot part of the age range
        end_age = end_age
        begin_age = window
        
        #Can be set to any interval, but the window shows the motion path that was 
        #calculated for the model
        interval = window
        
        #Create the time range for the motion path
        times = np.arange(begin_age,end_age,interval)
        
        
        iterations_array = np.arange(0,iterations,1)
        iteration_check = [25,50,75,100,125,150,175]
        
        #We first load each rotation file and compute all rotation poles, this
        #way the rotation file has to be read only once per iteration    
        for iteration in iterations_array:
            if iteration in iteration_check:
                print('iteration: ',iteration)
                
            dataframe_csv = pd.read_csv(output_dir+'/data/' + rot_file_folder + '.csv', sep=';')
            
            dataframe_iter = dataframe_csv.loc[dataframe_csv["iteration"]==iteration]
            
            age_range = []
            
            for time in dataframe_iter["age"]:
                if time < (end_age-(window/2)):
                    age_range.append(time)
                
            rotation_file_name = rot_file_inputdir +'/'+ rot_file_folder + '/optimised_rotation_model_' + rot_file_folder + str(iteration) + 'test.rot'
        
            rotation_features = pygplates.FeatureCollection(rotation_file_name)
        
            rotation_model = pygplates.RotationModel(rotation_file_name)
            
            for time in times:
                
                #Get the actual time of the iteration
                diff = 10
                actual_time = 0
                for time2 in age_range:
                    if abs(time2-time) <diff:
                        diff = abs(time2-time)
                        actual_time = time2
                        
                #Get the rotation of the plate from 0 Ma to time
                rotation = rotation_model.get_rotation(actual_time, plate_id)
                
                lat_pole0 = list(lat_pole0)
                lon_pole0 = list(lon_pole0)
                
                lat, lon, angle = rotation.get_lat_lon_euler_pole_and_angle_degrees()
                
                lat_c,lon_c = geoTools.checkLatLon(lat, lon)
                
                if test == True:
                    angle = -angle
                else:
                    angle = angle
                    
                euler_pole = [lat_c,lon_c,angle]
                euler_pole = list([euler_pole[0],euler_pole[1],euler_pole[2]])
                #print(euler_pole)
                
                pole = pmag.pt_rot(euler_pole,lat_pole0,lon_pole0)
                
                
                
                lat_pole = np.float64(pole[0][0])
                lon_pole = np.float64(pole[1][0])
                
                pole_dictionary.append({'lat': lat_pole, 'lon': lon_pole, 
                                        'time': time, 'iteration': iteration, 
                                        'model': label, 'plateid': plate_id,
                                        'actual_time': actual_time,
                                        'lat_pole': lat_pole0[0],'lon_pole':lon_pole0[0]})
                
                
                
                
                
                #pole_df = pole_df.append(new_row,ignore_index=True)
print(pole_dictionary)
new_df = pd.DataFrame(pole_dictionary)
print(new_df)

#pole_df_new = pd.concat([pole_df,new_df],ignore_index=True)

new_df.to_csv(output_dir+'/data/Motion_paths.csv', mode='a', header=False)

#%% Make a csv of mean poles describing absolute motion of Africa

#Coordinates of the point, default is set in Central Africa
lat_pole0 = list([-2.31])
lon_pole0 = list([19.15])

#Plate id of the plate where the point is located
plate_id = 701

pole_df_new = pd.read_csv(output_dir+'/data/Motion_paths.csv')

#Relative plate model
original_rot_file = os.path.join('data') + '/1000_0_rotfile_Merdith_et_al.rot'

#Remove any existing rotations with 701 as the moving plate
original_rotations = pygplates.FeatureCollection(original_rot_file)

filtered_features = []

for feature in original_rotations:
    total_reconstruction_pole = feature.get_total_reconstruction_pole()
    if total_reconstruction_pole:
        fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
        if moving_plate_id == 701:
            continue
        else:
            filtered_features.append(feature)
                

filtered__rotation_feature_collection = pygplates.FeatureCollection(filtered_features)


for j,age_error in enumerate(age_error_rot_files):
    
    rot_file_folder, colour, label, window, end_age, iterations,marker = age_error

    print(label)
    
    #Can be used to plot part of the age range
    end_age = 1000
    begin_age = window
    
    #Can be set to any interval, but the window shows the motion path that was 
    #calculated for the model
    interval = window
    
    #Create the time range for the motion path
    times = np.arange(begin_age,end_age,interval)
    
    mean_lats = []
    mean_lons = []
    
    #colormap to color the ellipses
    cmap2 = plt.get_cmap('viridis', int(end_age/interval))
    
    rotations = []
    
    rotations_feature_collection = []
    
    rotation_time_samples_701_rel_000 = []
    
    for n,time in enumerate(times):
        print(time,'Ma')
        
        df_lat_pole = pole_df_new.loc[pole_df_new["lat_pole"]==lat_pole0[0]]
        
        df_model = df_lat_pole.loc[df_lat_pole["model"]==label]
        
        
        dataframe_time = df_model.loc[df_model["time"]==time]
        
        dataframe_plateid = dataframe_time.loc[dataframe_time['plateid']==plate_id]
        
        lats = np.array(dataframe_plateid['lat'])
        lons = np.array(dataframe_plateid['lon'])
        
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
        
        finite_rotation = pygplates.FiniteRotation.create_great_circle_point_rotation((lat_pole0[0],lon_pole0[0]), (mean_lat,mean_lon))
        
        rotations_feature_collection.append(finite_rotation)
        
        rotation_time_samples_701_rel_000.append(
            pygplates.GpmlTimeSample(
                    pygplates.GpmlFiniteRotation(finite_rotation),
                    time,
                    label))
        
        pole_latitude, pole_longitude, angle_degrees = finite_rotation.get_lat_lon_euler_pole_and_angle_degrees()
        
        rotations.append({'Age': time, 'Latitude': np.round(pole_latitude,2), 
                          'Longitude': np.round(pole_longitude,2), 'Angle': np.round(angle_degrees,2)})
        
    df = pd.DataFrame(rotations)
    
    df.to_csv(output_dir+'/data/Poles_'+label+'.csv', sep = ';', index=False)
    
    optimised_rotation_feature = pygplates.Feature.create_total_reconstruction_sequence(
        0,
        701,
        pygplates.GpmlIrregularSampling(rotation_time_samples_701_rel_000))
    
    optimised_rotation_feature_collection = pygplates.FeatureCollection(optimised_rotation_feature)
    
    optimised_rotation_feature_collection.add(filtered__rotation_feature_collection)
    
    # filtered_features = []
    
    # for feature in optimised_rotation_feature_collection:
    #     total_reconstruction_pole = feature.get_total_reconstruction_pole()
    #     if total_reconstruction_pole:
    #         fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
    #         if fixed_plate_id == 70 and moving_plate_id == 701:
    #             continue
    #         else:
    #             filtered_features.append(feature)
                    
    
    # filtered__rotation_feature_collection = pygplates.FeatureCollection(filtered_features)    

    optimised_rotation_feature_collection.write(output_dir+'/data/'+label+'.rot')     
    
    
    # for rotation_feature_collection in optimised_rotation_feature_collection:
    
    #     rotation_feature_index = 0
        
    #     while rotation_feature_index < len(rotation_feature_collection):
    #         rotation_feature = rotation_feature_collection[rotation_feature_index]
    #         total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
    #         if total_reconstruction_pole:
    #             fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
                
    #             # Change fixed plate ID 000 to 005 (unless it's the 005-000 plate pair).
    #             # We want everything that references 000 to now reference 005.
    #             # Later we'll add a 005-000 sequence to store optimised rotation adjustments.
    #             if fixed_plate_id == 70:
    #                 if moving_plate_id == 701:
    #                     # Exclude 005-000 rotation features (we'll add a new one later).
    #                     del rotation_feature_collection[rotation_feature_index]
    #                     rotation_feature_index -= 1
            
    #         rotation_feature_index += 1
    
    
    # optimised_rotation_feature_collection.write(output_dir+'/data/'+label+'.rot')
        


#%% Make the motion path figure and uncertainty size over time figure

#If True, all points are plotted per timestep. This is useful when plotting a 
#specific timestep to see the spread of points, but becomes very cluttered
#when all timesteps are plotted
scatter = False

#Coordinates of the point, default is set in Central Africa
lat_pole0 = list([-2.31])
lon_pole0 = list([19.15])

pole_df_new = pd.read_csv(output_dir+'/data/Motion_paths.csv')

major_axis_length = []

for j,age_error in enumerate(age_error_rot_files):
    
    rot_file_folder, colour, label, window, end_age, iterations,marker = age_error
    
    
    
    print(label)
    
    #Can be used to plot part of the age range
    end_age = 350
    begin_age = window
    
    #Can be set to any interval, but the window shows the motion path that was 
    #calculated for the model
    interval = window
    
    #Set extent can be turned on or off
    set_extent = True
    
    #Create the time range for the motion path
    times = np.arange(begin_age,end_age,interval)
    
    
    #### Set up the figure ####
    
    if end_age >500:
        fig = plt.figure(figsize=(25,32),dpi=300)
    else:
        fig = plt.figure(figsize=(25,16),facecolor='white',dpi=300)
    
    #Set the lat and lon grid to be plotted on the map
    lat_grid=[90,80,60,40,20,0,-20.,-40.,-60.,-80.,-90.]
    lon_grid=[-180., -150., -120.,  -90.,  -60.,  -30.,    0.,   30.,   60., 90.,  120.,  150.,  180.]
    
    map_projection = ccrs.Orthographic(
        central_longitude=lon_pole0[0], central_latitude=lat_pole0[0])
    
    map_axis = fig.add_subplot(1,1,1, projection=map_projection)
    
    

    map_axis.set_global()


    map_axis.add_feature(cartopy.feature.LAND, zorder=0,
                        facecolor='burlywood', edgecolor='black')

    map_axis.add_feature(cartopy.feature.OCEAN, zorder=0, facecolor='lightcyan',alpha=0.6)

    gl = map_axis.gridlines(xlocs=lon_grid, ylocs=lat_grid, linewidth=1,
                      color='black', linestyle='solid',draw_labels=True)
    
    gl.xlabel_style = {'size': 40}
    gl.ylabel_style = {'size': 40}
    
    
    cmap = plt.get_cmap('viridis')
    
    # #Add in the colorbar
    # ticks = np.arange(begin_age-interval,end_age+100,100)
    # sm = plt.cm.ScalarMappable(
    #             cmap='viridis', norm=BoundaryNorm(range(begin_age-interval,end_age+10,10),cmap.N))
    # cb = plt.colorbar(sm, orientation = 'horizontal', ticks = ticks, shrink=0.6, pad=0.07,ax=map_axis) 
    # cb.set_label(label = 'Age (Ma)',size=45,labelpad=30)
    # cb.ax.set_xticklabels(ticks)
    # cb.ax.tick_params(labelsize=40,pad=10)
    

    mean_lats = []
    mean_lons = []
    
    #colormap to color the ellipses
    cmap2 = plt.get_cmap('viridis', int(end_age/interval))
    
    l1s = []
    for n,time in enumerate(times):
        print(time,'Ma')
        
        df_lat_pole = pole_df_new.loc[pole_df_new["lat_pole"]==lat_pole0[0]]
        
        df_model = df_lat_pole.loc[df_lat_pole["model"]==label]
        
        
        dataframe_time = df_model.loc[df_model["time"]==time]
        
        dataframe_plateid = dataframe_time.loc[dataframe_time['plateid']==plate_id]
        
        lats = np.array(dataframe_plateid['lat'])
        lons = np.array(dataframe_plateid['lon'])
        
        
        

        #Calculate the kent distribution to find the mean pole per timestep and
        #construct a uncertainty ellipse encompassing 95% of the datapoints
        kent_distribution_95 = ipmag.kent_distribution_95(dec=lons,inc=lats) 
        
        [l1,        #Length of the major axis of the ellipse
         l2,        #Length of the minor axis of the ellipse
         bearing,   #Azimuth of the major axis relative to N
         mean_lon,  #Mean longitude
         mean_lat,  #Mean latitude
         ] = ellipse_from_kent(kent_distribution_95)
        
        l1s.append(l1)
        
        
        
        
        mean_lats.append(mean_lat)
        mean_lons.append(mean_lon)
           
        #Determine zorder relative to time, this way the younger datapoints will
        #plot on top of the older datapoints
        zorder = len(times)-n
        
        #Plot fill of ellipse
        plot_pole_ellipse(map_axis,kent_distribution_95, 
                  color=cmap2(n), edgecolor=cmap2(n), marker='s', 
                  markersize=20, label='', lw=1, lower=False,filled=True, 
                  facecolor=cmap2(n), alpha=0.1, zorder=zorder)
        
        #Plot outline of ellipse
        plot_pole_ellipse(map_axis,kent_distribution_95, 
                  color=cmap2(n), edgecolor=cmap2(n), marker='s', 
                  markersize=20, label='', alpha=1, lw=1, lower=False,
                  zorder=zorder-1,filled=False)
        
        #Plot the mean poles per timestep
        map_axis.scatter(mean_lon,mean_lat, transform=ccrs.Geodetic(), 
                         c=cmap2(n),s=200,edgecolor=colour, marker=marker, 
                         zorder=zorder+120)
        
        #If scatter is True, all poles for certain timesteps will be plotted
        if scatter:
            if  time == 550:
                map_axis.scatter(lons,lats, transform=ccrs.Geodetic(), c=cmap2(n),
                                 s=5, zorder=2,alpha=0.5)
    
    major_axis_length.append(l1s)
    
    #Insert the seed point at 0 Ma so the motion path starts at 0 Ma
    if begin_age-10 == 0:
        mean_lons.insert(0,lon_pole0[0])
        mean_lats.insert(0,lat_pole0[0])
        times_plot = np.arange(0,end_age+interval,interval)
        
    else:
        times_plot = times
    
    #Plot the motion path as a line 
    map_axis.plot(mean_lons,mean_lats,transform=ccrs.Geodetic(),color = colour, 
                  lw=10,label=label,marker=marker,alpha=0.7,zorder=110)
    
    if set_extent:
        if end_age <500:
            map_axis.set_extent([-40,25,5,-20])
            
        else: 
            map_axis.set_extent([-40,40,60,-60])
    
    plt.tight_layout()
    
    if scatter:
        plt.savefig(output_dir+'/figures/' +'Motion_path_of_'+str(plate_id)+
                    '_in_the_'+label+'_'+str(begin_age)+'-'+str(end_age)+'Ma_with_scatter.png', dpi=300) 
    else:
        plt.savefig(output_dir+'/figures/' +'Motion_path_of_'+str(plate_id)+
                    '_in_the_'+label+'_'+str(begin_age)+'-'+str(end_age)+'Ma.png', dpi=300) 
        
    
 
#### Make the figure showing the ellipse size over time ####

fig,ax1 = plt.subplots(1,1,figsize=(40,10), dpi=100)
        
for j,age_error in enumerate(age_error_rot_files):
    
    rot_file_folder, colour, label, window, end_age, iterations,marker = age_error
    
    
    l1s = major_axis_length[j]
    
    times = np.arange(window,end_age,window)
    
    ax1.plot(times,l1s,c=colour, linewidth=6,alpha=1)
    
    ax1.scatter(times,l1s,c=colour,alpha=1)
    
    ax1.set_ylabel('Long axis in km',fontsize=40)
    ax1.set_xlabel('Time in Ma',fontsize=40,labelpad=20)
    ax1.set_ylim(0,3000)
    ax1.set_xlim(0,1000)
    ax1.grid(visible=True, color='darkgrey', linewidth=2)
    ax1.set_xticks(ticks=np.arange(0,end_age,100))
    ax1.tick_params(labelsize=30)
    ax1.set_title('Size of uncertainty ellipse', fontsize=60,fontweight='bold')
    
plt.savefig(output_dir+'/figures/' +'Size_of_major_axis_uncertainty_ellipse_in_the_'+label+'_'+str(begin_age)+'-'+str(end_age)+'Ma.png', dpi=300) 
  
        
#%%############################################################################
#Compute the net true polar wander angle over time and plot the TPWP as well 
#as a APWP of the Vaes 2023 GAPWAP and a synthetic APWP of the mantle reference
#frame computed here
###############################################################################
"""
This code has been minimally adapted from the code published in Vaes and
van Hinsbergen, 2024


"""

#Load Vaes 2023 GAPWAP
t_max = 320 # define maximum age 
dt = 10 # define time step of polar wander paths

df = Vaes_APWP_df

# Store key data in lists
lats_p,lons_p,ages=df['plat'],df['plon'],df['age']
N_ages=len(lats_p)


#%% Compute the TPWP path and net angle   
#Check if the data file already exists otherwise we create it here
file_check = os.path.exists(output_dir+'/data/TPW.csv') 

if file_check == False:
    print('File does not exist, making one now...')
    pole_df = pd.DataFrame(columns = ['lat', 'lon', 'age','iteration','model','net_TPW','actual_time'])
    pole_df.to_csv(output_dir+'/data/TPW.csv')
    
if file_check:
    print('File exists, loading now...')
    pole_df = pd.read_csv(output_dir+'/data/TPW.csv')
    

    

    

pole_dictionary = []

for j,age_error in enumerate(age_error_rot_files):
    
    
    rot_file_folder, colour, label, window, end_age, iterations,marker = age_error
    
    
    
    print(label)
    
    #Check if the input data is missing from the file
    missing_data = False
        
    if label in list(pole_df['model']):
        pole_df = pole_df.loc[pole_df["model"]==label]
        print('model')
    else:
        missing_data = True
    
    if missing_data == False:
        print('All data was computed before')
        
    elif missing_data:
        print('There is data missing! Calculating now...')   
        
        #Can be used to plot part of the age range
        end_age = end_age
        begin_age = window
        
        #Can be set to any interval, but the window shows the motion path that was 
        #calculated for the model
        interval = window
        
        #Create the time range for the motion path
        times = np.arange(begin_age,end_age,interval)
        
        
        iterations_array = np.arange(0,iterations,1)
        iteration_check = [25,50,75,100,125,150,175]
        
        #We first load each rotation file and compute all rotation poles, this
        #way the rotation file has to be read only once per iteration    
        for iteration in iterations_array:
            if iteration in iteration_check:
                print('iteration: ',iteration)
                
            dataframe_csv = pd.read_csv(output_dir+'/data/' + rot_file_folder + '.csv', sep=';')
            
            dataframe_iter = dataframe_csv.loc[dataframe_csv["iteration"]==iteration]
            
            age_range = []
            
            for time in dataframe_iter["age"]:
                if time < (end_age-(window/2)):
                    age_range.append(time)
                
            rotation_file_name = rot_file_inputdir +'/'+ rot_file_folder + '/optimised_rotation_model_' + rot_file_folder + str(iteration) + 'test.rot'
        
            rotation_features = pygplates.FeatureCollection(rotation_file_name)
        
            rotation_model = pygplates.RotationModel(rotation_file_name)
            
            rel_plate = 701
            fixed_plate = 0

            dic = []
            
            data = [[701,0,0,0,0,0],[701,330,0,0,0,0]]
            
            plate_IDs = [701,701]
            
            for time in times:
                
                #Get the actual time of the iteration
                diff = 10
                actual_time = 0
                for time2 in age_range:
                    if abs(time2-time) <diff:
                        diff = abs(time2-time)
                        actual_time = time2
                        
                #Get the rotation of the plate from 0 Ma to time
                rotation = rotation_model.get_rotation(actual_time, plate_id)
                
                lat, lon, ang = rotation.get_lat_lon_euler_pole_and_angle_degrees()
                
                ang = -ang
                
                data.append([1000,time,lat,lon,ang,701])
                
            for data_list in data:
                data_list = np.array(data_list)
                
            data = np.array(data)
            
            plate_IDs=data[:,0]     # stores plate IDs in array
            ages=data[:,1]          # stores ages of rotation poles in array
            file_length=len(plate_IDs)  # computes length of rotation file
            
            #--------------------
            # SPECIFY PLATE AND AGES
            
            plate_ID = 701 # plate ID of South Africa
            fixed_plate_ID = 1000 # plate ID of mantle
            
            init_plate_ID=plate_ID
            plate_ID_in_file=False
            for i in range(file_length):    # checks if plate ID is included in rotation file
                if (int(plate_IDs[i])==plate_ID):
                    
                    plate_ID_in_file=True
                    break
            if plate_ID_in_file==False:
                print('Error: Plate ID is not in file')
            
            if (fixed_plate_ID==plate_ID):
                print('Error: plate IDs are equal')
            fixed_plate_ID_in_file=False
            for i in range(file_length):    # checks if plate ID is included in rotation file
                if (plate_IDs[i]==fixed_plate_ID or fixed_plate_ID==0):
                    fixed_plate_ID_in_file=True
                    break
            if fixed_plate_ID_in_file==False:
                print('Error: Fixed plate ID is not in file')
                          
            # CREATE LIST OF AGES
            t = []
            
            for k in range(len(df['age'])):
                if df['age'][k]<=t_max:
                    t.append(df['age'][k])
                else:
                    break
            #print('Ages for TPWP:', t)
            
            # COMPUTE TOTAL RECONSTRUCTION POLES AT GIVEN TIMES
            tr_poles = []
            for i in range(len(t)):         # loops over selected times
                if (t[i]==0):           
                    print(0,',',0,',',0,',',t[i])    # prints pole at t=0 Ma
                    tr_poles.append([0,0,0])
                else:
                    plate_pole=rf.rot2frame(data,plate_ID,t[i])                 # computes pole of selected plate rel. to reference frame
                    fixed_plate_pole=rf.rot2frame(data,fixed_plate_ID,t[i])     # computes pole of fixed plate rel. to reference frame
                    fixed_plate_pole[2]=-1*fixed_plate_pole[2]          # changes sign of rotation angle
                    POLE=rf.addpoles(plate_pole,fixed_plate_pole)          # computes total reconstruction pole
                    tr_poles.append(POLE)
                    
            rlats,rlons=[],[]
            for i in range(len(t)):
                rlon,rlat=rf.pt_rot(tr_poles[i],[lats_p[i]],[lons_p[i]])
                rlats.append(rlat[0])
                rlons.append(rlon[0])
            
            # change sign of pole latitudes
            rlons = [rlons[i]-180 for i in range(len(tr_poles))]
            rlats = [rlats[i]*-1 for i in range(len(tr_poles))]
            
            rot_APWP=pd.DataFrame(list(zip(rlons,rlats,t)),columns=['plon','plat','age'])
            
            TPW_angles = []
            for i in range(len(rot_APWP)):
                TPW_angle = pmag.angle([rot_APWP['plon'][i],rot_APWP['plat'][i]],[0,90])
                TPW_angles.append(TPW_angle[0])
            
            for i in range(len(rlons)):
                pole_dictionary.append({'lat': rlats[i], 'lon': rlons[i], 'age': t[i],
                                   'iteration': iteration, 'model': label,
                                   'actual_time': actual_time, 'net_TPW': TPW_angles[i]})
                
new_df = pd.DataFrame(pole_dictionary)
print(new_df)

pole_df_new = pd.concat([pole_df,new_df],ignore_index=True)

pole_df_new.to_csv(output_dir+'/data/TPW.csv', mode='a', header=False)


#%%Make the figure
"""
Makes TPWPs, a graph with the net angle of TPW over time and an APWP of Vaes
et al 2023 and the mantle reference frame. Make sure the motion paths datafile 
contains the south pole coordinates for the APWP.

"""

pole_df_new = pd.read_csv(output_dir+'/data/TPW.csv')

pole_df_mp = pd.read_csv(output_dir+'/data/Motion_paths.csv')

for j,age_error in enumerate(age_error_rot_files):
    
    rot_file_folder, colour, label, window, end_age, iterations,marker = age_error 
    
    ### Plot the TPWP ###

    plt.rcParams['axes.facecolor']='white'
    #cmap = plt.cm.viridis
    
    lat_grid=[90,80,60,40,20,0,-20.,-40.,-60.,-80.,-90.]
    lon_grid=[-180., -150., -120.,  -90.,  -60.,  -30.,    0.,   30.,   60., 90.,  120.,  150.,  180.]

    fig = plt.figure(figsize=(40,32),facecolor='white',dpi=300)
    map_projection = ccrs.Orthographic(
        central_longitude=0, central_latitude=90)
    
    map_axis = fig.add_subplot(1,1,1, projection=map_projection)

    map_axis.set_global()

    
    map_axis.add_feature(cartopy.feature.LAND, zorder=0,
                        facecolor='burlywood', edgecolor='black')

    map_axis.add_feature(cartopy.feature.OCEAN, zorder=0, facecolor='lightcyan',alpha=0.6)

    map_axis.gridlines(xlocs=lon_grid, ylocs=lat_grid, linewidth=1,
                      color='black', linestyle='solid',draw_labels=True)
    
    cmap = plt.get_cmap('viridis')
    
    #Add in the colorbar
    ticks = [0,40,80,120,160,200,240,280,320]
    sm = plt.cm.ScalarMappable(
                cmap='viridis', norm=BoundaryNorm(range(0,330,10),cmap.N))
    cb = plt.colorbar(sm, orientation = 'horizontal', ticks = ticks, shrink=0.6, pad=0.06,ax=map_axis) 
    cb.set_label(label = 'Age (Ma)',size=40,labelpad=20)
    cb.ax.set_xticklabels(ticks)
    cb.ax.tick_params(labelsize=36,pad=10)      
    
    #Can be used to plot part of the age range
    end_age = t_max
    begin_age = 10
   
    interval = 10
    
    #Create the time range for the motion path
    times = np.arange(begin_age,end_age+interval,interval)
    
    mean_lons = []
    mean_lats = []
    
    cmap2 = plt.get_cmap('viridis', 32)
    
    min_ang_diff = []
    max_ang_diff = []
    APWP_diff = []
            
    for n,time in enumerate(times):
        
        print(time,'Ma')
        
        
        df_model = pole_df_new.loc[pole_df_new["model"]==label]
        
        
        dataframe_time = df_model.loc[df_model["age"]==time]
        
        # lons = list(dataframe_time.loc[dataframe_time['age']==time]['lon'])
        # lats = list(dataframe_time.loc[dataframe_time['age']==time]['lat'])
        lons = list(df_model.loc[df_model['age']==time]['lon'])
        lats = list(df_model.loc[df_model['age']==time]['lat'])
        
        net_TPW = np.array(df_model.loc[df_model['age']==time]['net_TPW'])
        min_ang_diff.append(np.min(net_TPW))
        max_ang_diff.append(np.max(net_TPW))
        
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
        
        #Calculate the mean TPW net angle
        TPW_angle = pmag.angle([mean_lon,mean_lat],[0,90])
        APWP_diff.append(TPW_angle[0])

        
        #Determine zorder relative to time, this way the younger datapoints will
        #plot on top of the older datapoints
        zorder = len(times)-n
        
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
        
    mean_lons.insert(0,0)
    mean_lats.insert(0,90)
    times_plot = np.arange(0,end_age+interval,interval)
    
    #Plot the APWP as a line 
    map_axis.plot(mean_lons,mean_lats,transform=ccrs.Geodetic(),color = colour, 
                  lw=10,label=label,marker=marker,alpha=0.7,zorder=110)
    
    plt.tight_layout()
    
    
    plt.savefig(output_dir+'/figures/' +'TPWP_in_the_'+label+'.png', dpi=300) 
        
    ## Plot TPW net angle ##
    
    
    fig = plt.figure(figsize=(25,11),dpi=150)
    
    plt.xlim(0,320)
    

     
    plt.grid(visible=True, color='darkgrey', linewidth=2)
    plt.scatter(times,APWP_diff, c=colour,zorder=4,edgecolors='k')
    plt.plot(times,APWP_diff,label=label, c=colour)
    plt.fill_between(times, min_ang_diff, max_ang_diff, color = 'teal',alpha=0.4)
    
    plt.plot(times,min_ang_diff,color=colour)
    plt.plot(times,max_ang_diff,color=colour)
    
    
    plt.ylabel('Net true polar wander in degrees', fontsize = 40,labelpad=20)
    plt.xlabel('Time in Ma', fontsize = 40,labelpad=20)
    
    plt.tick_params(labelsize=40)

    plt.ylim(0,42)
    
    #plt.title('Angular difference with Pmag frame', fontsize = 50)
    
    plt.savefig(output_dir+'/figures/' +'TPW_angle_in_the_'+label+'.png', dpi=300) 
    
###### Plot both APWPs for comparison ##########################
    
    lat_grid=[90,80,60,40,20,0,-20.,-40.,-60.,-80.,-90.]
    lon_grid=[-180., -150., -120.,  -90.,  -60.,  -30.,    0.,   30.,   60., 90.,  120.,  150.,  180.]

    fig = plt.figure(figsize=(40,32),facecolor='white',dpi=300)
    map_projection = ccrs.Orthographic(
        central_longitude=60, central_latitude=-50)

    cmap = plt.get_cmap('viridis')


    map_axis = fig.add_subplot(1,1,1, projection=map_projection)

    map_axis.set_global()


    map_axis.add_feature(cartopy.feature.LAND, zorder=0,
                        facecolor='burlywood', edgecolor='black')

    map_axis.add_feature(cartopy.feature.OCEAN, zorder=0, facecolor='lightcyan',alpha=0.6)

    map_axis.gridlines(xlocs=lon_grid, ylocs=lat_grid, linewidth=1,
                      color='black', linestyle='solid',draw_labels=True)
    
    #Add in the colorbar
    ticks = [0,40,80,120,160,200,240,280,320]
    sm = plt.cm.ScalarMappable(
                cmap='viridis', norm=BoundaryNorm(range(0,330,10),cmap.N))
    cb = plt.colorbar(sm, orientation = 'horizontal', ticks = ticks, shrink=0.6, pad=0.06,ax=map_axis) 
    cb.set_label(label = 'Age (Ma)',size=40,labelpad=20)
    cb.ax.set_xticklabels(ticks)
    cb.ax.tick_params(labelsize=36,pad=10) 
    
    #map_axis.set_extent([-30,100,-90,-30])
     
    plt.title('APWP for '+label, fontsize=50)
    
    ### Add the Vaes 2023 APWP ###
        
    APWP = df

    #APWP.loc[APWP['plon'] > 180, 'plon'] = APWP['plon'] - 360
    for i in range(len(APWP['plon'])):
        if np.float64(APWP['plon'][i]) >180:
            np.float64(APWP['plon'][i]) ==  np.float64(APWP['plon'][i])-360
             

    min_age,max_age = 0,320
    timestep = 10
    times_vaes = np.arange(10, max_age+timestep, timestep)

    cmap2 = plt.get_cmap('viridis', 33)
    # PLOT NEW GLOBAL APWP
    #ipmag.plot_poles_colorbar(map_axis, APWP.plon, APWP.plat, APWP.V95, APWP.age, min_age, max_age, markersize=150,edgecolor='k',filled_pole=True,fill_alpha=0.4,outline=True,colorbar=False)
    for i in range(len(times_vaes)):
        colour2 =cmap2(i)
        ipmag.plot_pole(map_axis,APWP.plon[i],APWP.plat[i],APWP.V95[i],markersize=50,marker='o',filled_pole=True,color=colour2,outline=True,fill_alpha=0.4, fill_color=colour2) # no colorbar
        
        if times_vaes[i] == 100 or times_vaes[i] == 250:
            map_axis.scatter(APWP.plon[i],APWP.plat[i], transform=ccrs.Geodetic(), c='r',s=200,edgecolor=colour, zorder=4)
            
    # plot connecting lines
    plt.plot(APWP.plon, APWP.plat, transform=ccrs.Geodetic(), color = 'k', linestyle = '--', lw=4,label = 'V23') # plot connecting lines

    ### Add the mantle reference frame ###
    
    df_model = pole_df_mp.loc[pole_df_mp["model"]==label]
    df_lat = df_model.loc[df_model["lat_pole"]==-90]
    
    mean_lons = []
    mean_lats = []
    
    for n,time in enumerate(times_vaes):
        
        df_time = df_lat.loc[df_lat['time']==time]
        
        lats = np.array(df_time['lat'])
        lons = np.array(df_time['lon'])
        
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
        
        if time == 100 or time == 250:
            map_axis.scatter(mean_lon,mean_lat, transform=ccrs.Geodetic(), c='r',s=200,edgecolor=colour, zorder=zorder+121,marker=marker)
    
        
    mean_lons.insert(0,0)
    mean_lats.insert(0,-90)
    times_plot = np.arange(0,end_age+interval,interval)
    
    #Plot the APWP as a line 
    map_axis.plot(mean_lons,mean_lats,transform=ccrs.Geodetic(),color = colour, 
                  lw=10,label=label,marker=marker,alpha=0.7,zorder=110)
    
    plt.tight_layout()
    
    
    plt.savefig(output_dir+'/figures/' +'APWP_in_the_'+label+'.png', dpi=300) 
        
    