# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 13:52:48 2024

@author: 6408885
"""

"""
WRITTEN BY STEFANIA WAGENAAR

This file uses the outcome of a TRM with age error and computes the net lithosphere
rotation and trench kinematics, and reconstructs the eruption locations of kimberlite
locations and LIP locations. This file can either be run in parallel or serially,
but serially takes really long.


"""

import os
import sys
import numpy as np
import itertools 
import os.path
import glob
import pandas as pd
from subduction_convergence_ptt import subduction_convergence as sc
import net_rotation
import geoTools
from pmagpy import pmag


#Select here to run in parallel or not
parallel = True

if parallel:
    from mpi4py import MPI
    
    # It seems that if one process/rank raises an exception then we need to manually
    # kill the other MPI processes according to:
    #
    #   https://groups.google.com/forum/#!topic/mpi4py/RovYzJ8qkbc
    #
    # ...otherwise MPI Finalize (in the process that raised exception) will block waiting for
    # the other processes to finish, but they're waiting for input (gather) from the rank=0 process
    # resulting in a deadlock.
    #
    # This code was obtained from:
    #
    #   https://groups.google.com/forum/#!topic/mpi4py/ktAZWIfx8zI
    #
    # ...and is the easiest way to do this if we don't care about properly cleaning up the processes.
    #
    _excepthook = sys.excepthook
    def excepthook(t,v,tb):
        _excepthook(t,v,tb)
        if (not MPI.Is_finalized()
            and MPI.Is_initialized()):
            MPI.COMM_WORLD.Abort(1)
    sys.excepthook = excepthook
    
    mpi_comm = MPI.COMM_WORLD
    mpi_size = mpi_comm.Get_size()
    mpi_rank = mpi_comm.Get_rank()


if parallel == False or mpi_rank == 0:

    #######################################################
    # Define Input files
    #######################################################
    
    age_error_rot_files = [
        ("dataframePV_1000", 'k', 'Continent_frame',10,1000,187, 'o'),
        ("dataframeTM_mean", 'r', 'Trench_frame',10,1000,99,'o'),
        ]
    
    rot_file_inputdir = os.path.join('optimisation')
    
    #Choose one of the frames to run
    #The slurm file needs a number of iterations so this needs to be edited as well
    
    rot_file_folder, colour, label, window, end_age, iterations,marker = age_error_rot_files[0]
    
    print(label)
    
    #Define the topology and continental polygons filenames
    #We will load these in the function, because gpml files cannot be distributed
    #among parallel processes
    
    topology_filenames = glob.glob(os.path.join('data', 'Topologies', '*.gpml'))
    
    continent_base_filenames = [
            'shapes_continents_Merdith_et_al.gpml'
        ]
    continent_filenames = [os.path.join('data', base_filename)
            for base_filename in continent_base_filenames]
    
    sys.stdout.flush()
    
    #Load the kimberlite file
    
    dataframe = pd.read_csv(os.path.join('data')+'/Kimberlites.csv',sep=';')
    dataframe = dataframe.convert_dtypes(convert_floating=True)
    dataframe = dataframe.reset_index()
    dataframe = dataframe.drop(dataframe[dataframe['Age'] == 0].index)
    dataframe = dataframe.reset_index()

    lons = []
    lats = []
    ages = []

    plate_ids = np.zeros(len(dataframe))
    
    numbers = []
    
    #The kimberlite file of Torsvik 2010 has no plate ids, only country names
    #so this section assigns every kimberlite with a plate id

    for i in range(len(dataframe)):
        lons.append(np.float64(dataframe['Longitude'][i]))
        lats.append(np.float64(dataframe['Latitude'][i]))
        ages.append(np.float64(dataframe['Age'][i]))
        numbers.append(np.float64(dataframe['No.'][i]))
        
        if 'South_Africa' in dataframe['Country'][i] or 'Mozambique' in dataframe['Country'][i] or 'Botswana' in dataframe['Country'][i] or 'Tanzania' in dataframe['Country'][i] or 'Namibia' in dataframe['Country'][i] or 'Congo' in dataframe['Country'][i] or 'Swaziland' in dataframe['Country'][i] or 'Lesotho' in dataframe['Country'][i]  or 'Angola' in dataframe['Country'][i]:
            plate_ids[i] = 701
            
        elif 'Australia' in dataframe['Country'][i]:
            plate_ids[i] = 801
            
        elif 'India' in dataframe['Country'][i]:
            plate_ids[i] = 501
            
        elif 'Antarctica' in dataframe['Country'][i]:
            plate_ids[i] = 802
            
        elif 'Greenland' in dataframe['Country'][i]:
            plate_ids[i] = 102
            
        elif 'USA' in dataframe['Country'][i] or 'Canada' in dataframe['Country'][i]:
            plate_ids[i] = 101
            
        elif 'Sweden' in dataframe['Country'][i] or 'Finland' in dataframe['Country'][i]:
            plate_ids[i] = 302
            
        elif 'Venezuela' in dataframe['Country'][i] or 'Brazil' in dataframe['Country'][i] or 'Bolivia' in dataframe['Country'][i]:
            plate_ids[i] = 201
            
        elif 'Guinea' in dataframe['Country'][i] or 'Mauritania' in dataframe['Country'][i] or 'Liberia' in dataframe['Country'][i] or 'Mali' in dataframe['Country'][i] or 'Gabon' in dataframe['Country'][i] or 'Ghana' in dataframe['Country'][i]:
            plate_ids[i] = 714
            
        elif 'Russia' in dataframe['Country'][i]:
            plate_ids[i] = 401
            
        elif 'China' in dataframe['Country'][i]:
            plate_ids[i] = 601
        
    #Collect data per kimberlite
    kimberlites = []

    for i in range(len(lons)):
        lst = [lons[i],lats[i],ages[i],plate_ids[i],numbers[i]]
        kimberlites.append(lst)
    
    print('kimberlite data loaded')
    
    sys.stdout.flush()
    
    dataframe = pd.read_csv(os.path.join('data')+'/LIPs.csv',sep=';')
    dataframe = dataframe.convert_dtypes(convert_floating=True)
    dataframe = dataframe.reset_index()
    dataframe = dataframe.drop(dataframe[dataframe['Age'] == 0].index)
    dataframe = dataframe.reset_index()


    lons = []
    lats = []
    ages = []

    plate_ids = np.zeros(len(dataframe))
    
    numbers = []

    for i in range(len(dataframe)):
        lons.append(np.float64(dataframe['Long.'][i]))
        lats.append(np.float64(dataframe['Lat.'][i]))
        ages.append(np.float64(dataframe['Age'][i]))
        numbers.append(np.float64(dataframe['No.'][i]))
        plate_ids[i] = np.float64(dataframe['ID'][i])
        
    #Collect data per LIP
    LIPs = []

    for i in range(len(lons)):
        lst = [lons[i],lats[i],ages[i],plate_ids[i],numbers[i]]
        LIPs.append(lst)    
        
    print('LIP data loaded')
    
    sys.stdout.flush()
    
    
if parallel:
    
    if mpi_rank == 0:
    
        
        
        iteration_list = list(np.arange(0,iterations,1))

        iteration_list = iteration_list[0:-1]
        
        iterations = []
    
        for iteration in iteration_list:
            iterations.append([iteration])
            
        print(iterations)
        sys.stdout.flush()
            
        constantStartingConditions = [kimberlites, LIPs, rot_file_inputdir,
                                      rot_file_folder,topology_filenames, end_age,
                                      window, label]
            
        
        
    else:
        iterations = None
        constantStartingConditions = None
        
    #mpi_comm.barrier()
    
    
        
    iteration_list = mpi_comm.scatter(iterations, root=0)


    constantStartingConditions = mpi_comm.bcast(constantStartingConditions, root=0)
    
    iterations = []
    iterations.append(iteration_list)
    
    startingConditions = []
    startingConditions.append(iterations)
    startingConditions.extend(constantStartingConditions)
    
    (iterations,kimberlites, LIPs, rot_file_inputdir,rot_file_folder,
     topology_filenames, end_age, window, label) = startingConditions
    
    
    
def run_calc(iteration,kimberlites, LIPs, rot_file_inputdir,rot_file_folder,
             topology_filenames, end_age, window, label, nr_interval = 5):
    
    
    #If running serially, comment this
    iteration = iteration[0]
    
    #First do the trench calculations
    
    import numpy as np
    import pygplates
    import pandas as pd
    
    if end_age>980:
        end_age=980
        end_age = 50
    
    times = np.arange(window,end_age,window)
    
    tm_dictionary = []
    
    topology_features = pygplates.FeatureCollection()
    for topology_filename in topology_filenames:
        # (omit files with the string "inactive" in the filepath)
        if "Inactive" not in topology_filename:
            topology_features.add( pygplates.FeatureCollection(topology_filename) )
        else:
            topology_filenames.remove(topology_filename)
            
    
    
    rotation_file_name = rot_file_inputdir +'/'+ rot_file_folder + '/optimised_rotation_model_' + rot_file_folder + str(iteration) + 'test.rot'

    rotation_model = pygplates.RotationModel(rotation_file_name)
    
    rotation_features = pygplates.FeatureCollection(rotation_file_name)
    
    threshold_sampling_distance_radians = 1
    
    for k in range(len(times)):
        time = times[k]

        interval = window
        
        if iteration == 1:
            print(time)
            sys.stdout.flush()
            
        subduction_data = sc(
            rotation_model,
            topology_features,
            threshold_sampling_distance_radians,
            float(time),
            anchor_plate_id=0,
            velocity_delta_time=interval,
            output_subducting_absolute_velocity_components = True)
        
        subduction_data = np.vstack(subduction_data)
            
        subduction_lon     = subduction_data[:,0]
        subduction_lat     = subduction_data[:,1]
        subduction_angle   = subduction_data[:,3]
        subduction_norm    = subduction_data[:,7]
        subduction_pid_sub = subduction_data[:,8]
        subduction_pid_over= subduction_data[:,9]
        #subduction_length  = np.radians(subduction_data[:,6])*gplately.tools.geocentric_radius(subduction_data[:,1])
        subduction_convergence = np.fabs(subduction_data[:,2]) * np.cos(np.radians(subduction_data[:,3]))
        subduction_migration   = np.fabs(subduction_data[:,4]) * np.cos(np.radians(subduction_data[:,5]))
        #subduction_plate_vel = subduction_data[:,10]
        tm_vel_par = np.abs(subduction_data[:,11])*10

        # remove entries that have "negative" subduction
        # this occurs when the subduction obliquity is greater than 90 degrees
        #subduction_convergence = np.clip(subduction_convergence, 0, 1e99)
        

        
        
        
        trench_vel = subduction_data[:,4]*10
        trench_obl = subduction_data[:,5]
        slab_vel = subduction_data[:,2]*10
        slab_obl = subduction_data[:,3]
        tlon = subduction_lon
        tlat = subduction_lat
        arc_length = subduction_data[:,6]
        subducting_arc_normal_azimuth = subduction_data[:,7]
        subduction_zone_plate_id = subduction_pid_over
        
        
        tm_vel_orth = np.abs(trench_vel) * -np.cos(np.radians(trench_obl)) 
        tm_vel_orth = -tm_vel_orth
        
        
        sc_vel_orth = subduction_convergence*10
        
        
        
        #Calculate mean advancing velocity and mean retreating vel
        temp_1 = []
        temp_2 = []
        max_rollback = 0
        max_advance = 0
        max_parallel = 0
        id_rollback = 0
        id_advance = 0
        id_parallel = 0
        
        fast_rollback = []
        fast_advance = []
        for i in range(len(tm_vel_orth)):
            if np.abs(tm_vel_par[i]) > max_parallel:
                max_parallel = np.abs(tm_vel_par[i])
                id_parallel = subduction_zone_plate_id[i]
            #print(tm_vel_orth[i])
            if tm_vel_orth[i] >0:
                temp_1.append(tm_vel_orth[i])
                if tm_vel_orth[i] > 100:
                    fast_rollback.append(tm_vel_orth[i])
                if tm_vel_orth[i] > max_rollback:
                    max_rollback = tm_vel_orth[i]
                    id_rollback = subduction_zone_plate_id[i]
            if tm_vel_orth[i] <0:
                temp_2.append(tm_vel_orth[i])
                if tm_vel_orth[i] < -30:
                    fast_advance.append(tm_vel_orth[i])
                if tm_vel_orth[i] < max_advance:
                    max_advance = tm_vel_orth[i]
                    id_advance = subduction_zone_plate_id[i]
                
        rollback_vel = np.mean(temp_1)
        # r_std = np.std(temp_1)
        nr_rollback = len(temp_1)
        advance_vel = np.mean(temp_2)
        # a_std = np.std(temp_2)
        nr_advance = len(temp_2)
        nr_total = len(tm_vel_orth)
        
        rollback_perc = 100*(nr_rollback/nr_total)
        
        parallel_vel = np.mean(np.abs(tm_vel_par))
                
        tm_dictionary.append({'iteration': iteration, 'age': time, 'model': label, 'rollback_mean': rollback_vel,
                               'advance_mean': advance_vel, 'parallel_mean': parallel_vel, 'rollback_max': max_rollback,
                               'advance_max': max_advance, 'parallel_max': max_parallel, 'rollback_max_id': id_rollback,
                               'advance_max_id': id_advance, 'parallel_max_id': id_parallel, 'rollback_perc': rollback_perc,
                               'fast_rollback': 100*(len(fast_rollback)/nr_total), 'fast_advance': 100*(len(fast_advance)/nr_total)})
    
    #Then the net rotation calculations
    
    times_nr = np.arange(nr_interval,end_age,nr_interval)
    
    nr_dictionary = []
    
    for n,time in enumerate(times_nr):
        rot = net_rotation.calculate_net_rotation_internal_gplates(
                    rotation_features,
                    topology_features,
                    reconstruction_time = time,
                    #velocity_method = VelocityMethod.T_TO_T_MINUS_DT,
                    #velocity_method = VelocityMethod.T_PLUS_DT_TO_T,
                    #velocity_method = VelocityMethod.T_PLUS_MINUS_HALF_DT,
                    velocity_delta_time = nr_interval
                    )
            
        lat, lon, ang = rot.get_lat_lon_euler_pole_and_angle_degrees()
    
        nr_dictionary.append({'iteration': iteration, 'age': time, 'model': label, 'nr': ang, 'interval': nr_interval})
        
    #Then the kimberlites and LIPs
        
    dict_kimb = []
    dict_lips = []
    
    kimb_check = [1,100,200,300,400,500,600,700,800,900]
    
    for kimb in kimberlites:
        k_lon, k_lat, age, plate_id, number = kimb
        
        if iteration == 1 and number in kimb_check:
        
            print(number)
            sys.stdout.flush()
        
        rotation = rotation_model.get_rotation(age,int(plate_id))
        
        lat, lon, angle = rotation.get_lat_lon_euler_pole_and_angle_degrees()
        
        lat_c,lon_c = geoTools.checkLatLon(lat, lon)
        
        euler_pole = [lat_c,lon_c,angle]
        euler_pole = list([euler_pole[0],euler_pole[1],euler_pole[2]])
        #print(euler_pole)
        
        pole = pmag.pt_rot(euler_pole,[k_lat],[k_lon])
        
        lat_pole = np.float64(pole[0])
        lon_pole = np.float64(pole[1])
        
        dict_kimb.append({'iteration': iteration, 'lat': lat_pole[0], 'lon': lon_pole[0],
                          'age': age, 'number': number, 'plateid': plate_id})
        
        
        
    for lip in LIPs:
        lip_lon, lip_lat, age, plate_id, number = lip
        
        rotation = rotation_model.get_rotation(age,int(plate_id))
        
        lat, lon, angle = rotation.get_lat_lon_euler_pole_and_angle_degrees()
        
        lat_c,lon_c = geoTools.checkLatLon(lat, lon)
        
        euler_pole = [lat_c,lon_c,angle]
        euler_pole = list([euler_pole[0],euler_pole[1],euler_pole[2]])
        #print(euler_pole)
        
        pole = pmag.pt_rot(euler_pole,[lip_lat],[lip_lon])
        
        lat_pole = np.float64(pole[0])
        lon_pole = np.float64(pole[1])
        
        dict_lips.append({'iteration': iteration, 'lat': lat_pole[0], 'lon': lon_pole[0],
                          'age': age, 'number': number, 'plateid': plate_id})
    
    return tm_dictionary, nr_dictionary, dict_kimb, dict_lips


if parallel:
    for iteration in iterations:
        tm_dictionary, nr_dictionary, dict_kimb, dict_lips = run_calc(iteration,
                                                            kimberlites, LIPs, 
                                                            rot_file_inputdir,
                                                            rot_file_folder,
                                                            topology_filenames,
                                                            end_age, window, 
                                                            label)
        
    tm_dictionary = mpi_comm.gather(tm_dictionary, root=0)
    nr_dictionary = mpi_comm.gather(nr_dictionary, root=0)
    dict_kimb = mpi_comm.gather(dict_kimb, root=0)
    dict_lips = mpi_comm.gather(dict_lips, root=0)

    if mpi_rank == 0:

        tm_dictionary = list(itertools.chain.from_iterable(tm_dictionary))
        nr_dictionary = list(itertools.chain.from_iterable(nr_dictionary))
        
        df_tm = pd.DataFrame(tm_dictionary)
        df_nr = pd.DataFrame(nr_dictionary)
        
        dict_kimb = list(itertools.chain.from_iterable(dict_kimb))
        dict_lips = list(itertools.chain.from_iterable(dict_lips))
        
        df_kimb = pd.DataFrame(dict_kimb)
        df_lips = pd.DataFrame(dict_lips)
        
        
        
        sys.stdout.flush()
        
        df_tm.to_csv(os.path.join('Output')+'/data/TM_data'+label+'.csv', sep = ';')
        df_nr.to_csv(os.path.join('Output')+'/data/NR_data'+label+'.csv', sep = ';')
        
        df_kimb.to_csv(os.path.join('Output')+'/data/Kimberlites_rot'+label+'.csv', sep = ';')
        df_lips.to_csv(os.path.join('Output')+'/data/LIPs_rot'+label+'.csv', sep = ';')
        
  
else:
    iterations = np.arange(0,iterations,1)
    
    full_dict_tm = []
    full_dict_nr = []
    
    full_dict_kimb = []
    full_dict_lips = []
    
    for iteration in iterations:
        tm_dictionary, nr_dictionary, dict_kimb, dict_lips = run_calc(iteration,
                                                            kimberlites, LIPs, 
                                                            rot_file_inputdir,
                                                            rot_file_folder,
                                                            topology_filenames,
                                                            end_age, window, 
                                                            label)
        
        for entry in tm_dictionary:
            full_dict_tm.append(entry)
            
        for entry in nr_dictionary:
            full_dict_nr.append(entry)
            
        for entry in dict_kimb:
            full_dict_kimb.append(entry)
            
        for entry in dict_lips:
            full_dict_lips.append(entry)
        
    df_tm = pd.DataFrame(full_dict_tm)
    df_nr = pd.DataFrame(full_dict_nr)   
    
    df_tm.to_csv(os.path.join('Output')+'/data/TM_data'+label+'.csv', sep = ';')
    df_nr.to_csv(os.path.join('Output')+'/data/NR_data'+label+'.csv', sep = ';')
    
    df_kimb = pd.DataFrame(full_dict_kimb)
    df_lips = pd.DataFrame(full_dict_lips)   
    
    df_kimb.to_csv(os.path.join('Output')+'/data/Kimberlites_rot'+label+'.csv', sep = ';')
    df_lips.to_csv(os.path.join('Output')+'/data/LIPs_rot'+label+'.csv', sep = ';')

