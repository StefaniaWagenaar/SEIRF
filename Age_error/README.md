# **Age error alterations**

The main changes done to accomodate age error into the TRM was to loop the main running file and create a new age range each iteration that is randomly drawn from a window around timesteps. 
Each of the supporting files that use the age range or intervals have been edited to use the current interval rather than the set interval. At the end of the main running file, the rotation file is opened
and the rotation poles of 701 and 101 relative to 000 are pulled and saved in a csv. This dataframe can then be used as input in the Make_rot_files.py file to make rotation files for each iteration. 
These rotation files and the csv file with the dataframe can then be copied into the SEIRF folder to compute all the SEIRF components. 
