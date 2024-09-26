# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 15:34:58 2024

@author: 6408885
"""

"""
This file uses the dataframe from the age error code and makes rotation files 
for each iteration

"""

import numpy as np
import pygplates
import pandas as pd
import os

from optimised_rotation_updater_age_error import OptimisedRotationUpdater

def get_reference_params(age):
    """
    Returns a 2-tuple containg reference plate ID and reference rotation filename (or None).
    
    If reference rotation filename is None then it means the no-net-rotation model should be used.
    """
    
    ref_rotation_file = '1000_0_rotfile_Merdith_et_al'
    
    if age <= 550:
        ref_rotation_plate_id = 701
    else:
        ref_rotation_plate_id = 101
   
    
    return ref_rotation_plate_id, ref_rotation_file

def error_rotation_file(datadir,original_rotation_filenames,start_age,end_age,age_range,model_name, rotations, rank):
    
    data_model = 'Global_1000-0_Model_2017'
    
    optimised_rotation_updater = OptimisedRotationUpdater(
            datadir,
            original_rotation_filenames,
            start_age,
            end_age,
            age_range,
            get_reference_params,
            data_model,
            model_name,
            rank)
    
    
    
    
    for n,age in enumerate(age_range):
        
        
        
        ref_rotation_plate_id, _ = optimised_rotation_updater.reference_params_function(age)
        
        optimised_rotation_updater.update_optimised_rotation(rotations[n], ref_rotation_plate_id,age)
        
        #optimised_rotation_updater.save_to_rotation_files()
        
        filename = optimised_rotation_updater.optimised_rotation_filename[:39] + str(rank) + optimised_rotation_updater.optimised_rotation_filename[39:]
        
    
    return optimised_rotation_updater.optimised_rotation_filename


age_error_csvs_1000Ma = [
    #("Age_error_NR_1000Ma_Results_age_error_iterations", 'green', 'NR', 10, 1000, 36,'*'),
    ("dataframeTM_mean", 'darkorange', 'Trench frame',10,1000,98,'o'),
    ("dataframePV_1000", 'magenta', 'Continent frame',10,1000,187, 'o'),
    #("dataframeNR", 'green', 'NR', 10, 320, iterations[8]+iterations[9]),
    ]

age_error_csv = age_error_csvs_1000Ma[0]

csv_name, colour, label, window, end_age, iterations,marker = age_error_csv

dataframe_csv = pd.read_csv(os.path.join('dataframes')+'/'+csv_name+'.csv', sep=';')

iterations_array = np.arange(0,iterations,1)

it_check = [10,25,50,75,100,125,150,175]

original_rotation_filenames = ['Global_1000-0_Model_2017/1000_0_rotfile_Merdith_et_al.rot']

datadir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data', '')

model_name = csv_name

#model_dir = 'Global_Model_WD_Internal_Release_2019_v3'
data_model = 'Global_1000-0_Model_2017'
model_dir = 'Global_1000-0_Model_2017'


    

for iteration in iterations_array:
    if iteration in it_check:
        print(iteration)
        
    dataframe_iter = dataframe_csv.loc[dataframe_csv["iteration"]==iteration]

    
    times = np.arange(window,end_age,window)
    
    rotations = []
    
    age_range = []

    
    for time in dataframe_iter["age"]:
        if time < (end_age-(window/2)):
            age_range.append(time)
    
    for time in times:
        dataframe_time = dataframe_iter.loc[(dataframe_iter["age"]>=(time-(window/2))) & (dataframe_iter["age"]<(time+(window/2)))]

 
        if np.float64(dataframe_time["age"]) <=550:

            rotation = pygplates.FiniteRotation((np.float64(dataframe_time["lat701"]),np.float64(dataframe_time["lon701"])),np.radians(np.float64(dataframe_time["angle701"])))

            rotations.append(rotation)
            
        elif np.float64(dataframe_time["age"]) >550:
        
        
            
            rotation = pygplates.FiniteRotation((np.float64(dataframe_time["lat101"]),np.float64(dataframe_time["lon101"])),np.radians(np.float64(dataframe_time["angle101"])))
 
            rotations.append(rotation)
        
    
    
    rotation_file_name =  error_rotation_file(datadir,original_rotation_filenames,end_age,0,age_range,model_name, rotations, iteration)

