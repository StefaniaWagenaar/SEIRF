# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 16:34:51 2024

@author: 6408885
"""

"""
WRITTEN BY STEFANIA WAGENAAR

This file loads all data needed for the visualisation and computation of the 
properties of the SEIRF. 

"""


import os.path
import pandas as pd
from matplotlib import image
import gplately

output_dir = os.path.join('Output')
age_error_inputdir = 'C:/Users/6408885/optAPM-master/optAPM-master/model_output/Age_errors/'
rot_file_inputdir = os.path.join('optimisation')


age_error_csvs_1000Ma = [
    ("dataframePV_1000", 'k', 'Continent_frame',10,1000,187, 'o'),
    #("dataframeTM_mean", 'r', 'Trench_frame',10,1000,99,'o'),
    ]


age_error_rot_files = [
    ("dataframePV_1000", 'k', 'Continent_frame',10,1000,187, 'o'),
    ("dataframeTM_mean", 'r', 'Trench_frame',10,1000,99,'o'),
    ]


#### Paleomagnetic frame of Vaes 2023

Vaes_APWP_df = pd.read_csv(os.path.join('data')+'/Vaes_int_APWP_10Ma.csv',sep=',')


#### Hotspot data and etopo file

etopo1 = gplately.Raster(data=image.imread(os.path.join('data')+'/color_etopo1_ice_low.tif'))
etopo1.lats = etopo1.lats[::-1]








