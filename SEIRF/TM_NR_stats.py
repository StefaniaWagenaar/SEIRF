# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 16:25:30 2024

@author: 6408885
"""

"""
WRITTEN BY STEFANIA WAGENAAR


This file loads the datasets produced with file [name] and plots figures
showing the computed net rotation, trench kinematics with 
the associated uncertainty. 



"""

from DataLoader import *
from Functions import *

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



#%% Trench kinematics with age error
##############################################################################
## Plot the trench kinematics with age error
##############################################################################                                                                       


def plot_subplot(ax, colors, arrays, title, xlabel, ylabel, times, window, dataframe, ylim_min,ylim_max, double = True):
    
    """
    Function that makes a subplot of the desired parameter. In the code computing
    the trench kinematics with age errors, there are many possible parameters that
    can be exported so this function gives some flexibility in what is displayed. 
    
    
    colors: list of color strings for each parameter
    arrays: list of strings with the column name of the dictionary
    double: if True plots each interval as a block, rather than a point
    
    """
    
    #Make the time array for the plot, needs to be doubled if double is True
    
    if double:
        times_plot = [0]
        
        for time in times:
            times_plot.append(time)
            times_plot.append(time)
            
        times_plot = times_plot[0:-1]
        
    else:
        times_plot = times
        
    
    #Plot each parameter that was given as input     
        
    for i in range(len(arrays)):
        array = arrays[i]
        color = colors[i]
        
        #Make empty lists to put the minimum, maximum, and mean values in
        
        array_min = []
        array_mean = []
        array_max = []
        
        for time in times:
            #Locate the parameter values for the time
            dataframe_time = dataframe.loc[(dataframe["age"]>=(time-(window/2))) & (dataframe["age"]<(time+(window/2)))]
            
            dataframe_array = np.abs(np.array(dataframe_time[array]))
            
            if double:
            
                array_min.append(np.min(dataframe_array))
                array_mean.append(np.mean(dataframe_array))
                array_max.append(np.max(dataframe_array))
                array_min.append(np.min(dataframe_array))
                array_mean.append(np.mean(dataframe_array))
                array_max.append(np.max(dataframe_array))
                
            else:
                array_min.append(np.min(dataframe_array))
                array_mean.append(np.mean(dataframe_array))
                array_max.append(np.max(dataframe_array))
            
        #Plot the mean, minimum, and maximum value per time step    
        lw_extremes = 2  
            
        ax.plot(times_plot,array_min,c=color,lw=lw_extremes)
        ax.plot(times_plot,array_max,c=color,lw=lw_extremes)
        ax.plot(times_plot,array_mean,c=color,lw=4)
        
        #Fill between the minimum and maximum value
        ax.fill_between(times_plot, array_min,array_max,color=color,alpha=0.3)
        
        ax.set_title(title, fontsize=80, fontweight='bold')
        ax.set_xlabel(xlabel, fontsize = 60,labelpad=20)
        ax.set_ylabel(ylabel, fontsize = 60,labelpad=20)
        ax.grid(visible=True, color='darkgrey', linewidth=2)
        ax.set_ylim(ylim_min,ylim_max)
        ax.set_xlim(0,end_age)
        ax.tick_params(labelsize=50)      

    return  

#Iterate over the different frames

for n,age_error_csv in enumerate(age_error_csvs_1000Ma):
    
    csv_name, colour, label, window, end_age, iterations,marker = age_error_csv
    
    #Input the maximum age for the plot
    end_age = 350 #1000
    

    TMdataframe = pd.read_csv(os.path.join('Output')+'/data/TM_data_'+label+'.csv', sep = ';')
    
    times = np.arange(window,end_age+window,window)
    
    fig = plt.figure(figsize=(30,40))  
    
    fig_rows, fig_columns = 3,1
    
    #plt.suptitle('Trench kinematics in the '+label,fontsize=60, y=1)

    ax = fig.add_subplot(fig_rows,fig_columns,1)
    
    colors = ['blue', 'navy']
    arrays = ['rollback_mean','rollback_max']
    title = 'Orthogonal trench retreat'
    xlabel = 'Time in Ma'
    ylabel = 'Velocity in mm/yr'
    ylim_min = 0
    ylim_max = 200
    
    
    plot_subplot(ax, colors, arrays, title, xlabel, ylabel, times, window, TMdataframe, ylim_min,ylim_max, double = False)
    
    ax = fig.add_subplot(fig_rows,fig_columns,2)
    
    colors = ['red','maroon']
    arrays = ['advance_mean','advance_max']
    title = 'Orthogonal trench advance'
    xlabel = 'Time in Ma'
    ylabel = 'Velocity in mm/yr'
    ylim_min = 0
    ylim_max = 200
    
    
    plot_subplot(ax, colors, arrays, title, xlabel, ylabel, times, window, TMdataframe, ylim_min,ylim_max, double = False)
    
    ax = fig.add_subplot(fig_rows,fig_columns,3)

    
    colors = ['darkorchid','indigo']
    arrays = ['parallel_mean','parallel_max']
    title = 'Parallel slab dragging'
    xlabel = 'Time in Ma'
    ylabel = 'Velocity in mm/yr'
    ylim_min = 0
    ylim_max = 200
    
    
    plot_subplot(ax, colors, arrays, title, xlabel, ylabel, times, window, TMdataframe, ylim_min,ylim_max, double = False)

    
    
    
    # ax = fig.add_subplot(fig_rows,fig_columns,3)
    
    # colors = ['blue']
    # arrays = ['rollback_perc']
    # title = 'Percentage of trench segments retreating'
    # xlabel = 'Time in Ma'
    # ylabel = 'Percent of trench segments'
    # ylim_min = 0
    # ylim_max = 100
    
    
    # plot_subplot(ax, colors, arrays, title, xlabel, ylabel, times, dataframe, ylim_min,ylim_max, double = False)
    
    plt.tight_layout(pad=0.3)
    
    plt.savefig(output_dir+'/figures/' +'Trench_kinematics_'+label+str(end_age)+'.png', dpi=150)

#%% Net lithospheric rotation with age error
##############################################################################
## Plot the net lithospheric rotation with age error
############################################################################## 

for n,age_error_csv in enumerate(age_error_csvs_1000Ma):
    
    csv_name, colour, label, window, end_age, iterations,marker = age_error_csv
    
    #Input the maximum age for the plot
    end_age = 350 #1000
    begin_age = 0
    
    #Color to plot uncertainty
    color_nr = 'steelblue'
    
    
    dataframe_nr = pd.read_csv(age_error_inputdir + csv_name +'_output_nr2.csv', sep=';')
    
    nr_i_min = []
    nr_i_max = []
    nr_i_mean = []
    
    #Set the window to the interval window set when calculating the net 
    #rotation values. We choose to go with 5 Ma here following Atkins and Coltice
    #2021.
    window = 5
    

    times_interval = np.arange(window,end_age+window,window)
    
    
    
    for time in times_interval:
        #Locate the parameter values for the time
        dataframe_time = dataframe_nr.loc[dataframe_nr["interval"]==time]
        
        #Use this one if you want a graph with net rotation per interval
        nr_array = np.array(dataframe_time['nr_i']) 
        
        #Use this one if you want a graph with net rotation per Ma
        #nr_array = np.array(dataframe_time['nr_my'])
        
        nr_i_min.append(np.min(nr_array))
        nr_i_mean.append(np.mean(nr_array))
        nr_i_max.append(np.max(nr_array))
        
        
        
    max_150 = [1.3,1.3]
    min_150 = [0,0]
    time_150 = [150,350]
    
    time_83 = [83,350]
    
    range_min = [0.03,0.03]
    range_max = [0.25,0.25]
    time_range = [0,350]
    time_range_83 = [83,350]
    time_range_150 = [150,350]
    
    
        
    fig = plt.figure(figsize=(40,22),dpi=150)
    
    plt.fill_between(time_range,range_min,range_max,color = 'grey', alpha=0.4)
    plt.fill_between(time_range_83,range_min,range_max,color = 'grey', alpha=0.4)
    plt.fill_between(time_range_150,range_min,range_max,color = 'grey', alpha=1)
        
    plt.fill_between(times_interval, nr_i_min,nr_i_max, color = color_nr,label='net rotation per interval',alpha=1)
    plt.plot(times_interval,nr_i_min,color = color_nr)
    plt.plot(times_interval,nr_i_max,color = color_nr)
    plt.plot(times_interval,nr_i_mean,color=colour)
    plt.scatter(times_interval,nr_i_mean,color=colour)
    
    plt.fill_between(time_150,min_150,max_150,color = 'w', alpha=0.4)
    plt.fill_between(time_83,min_150,max_150,color = 'w', alpha=0.4)
    
    
        
    plt.ylim(0,1.05)
    plt.xlim(begin_age,end_age)

    plt.xlabel('Time in Ma', fontsize = 60,labelpad=20)
    plt.ylabel('Net rotation in degrees/Ma', fontsize = 60,labelpad=20)

    plt.tick_params(labelsize=50)

    plt.grid(visible=True, color='darkgrey', linewidth=2)

    plt.title('Net lithospheric rotation predicted by the '+label, fontsize = 70, fontweight='bold')

    plt.savefig(output_dir+'/figures/'+'Net_rotation_' +label+'.png', dpi=150)


