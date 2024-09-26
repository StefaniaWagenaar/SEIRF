# **Towards a Solid Earth Integrated Reference Frame (SEIRF)**

### **by S.D.M. Wagenaar, B. Vaes & D.J.J. van Hinsbergen**

This manuscript has been submitted to *Journal of Geophysical Research* in 2024:

-------
Repo for the data files and Python codes used to compute all components of the SEIRF. 

This repo contains the following sections:
- The Age_error folder contains Python files with the necessary alterations to compute the Tectonic Rules Model (TRM) as described in Tetley et al., 2019 with an age error. The files have the same names of the files they replace in the Earthbyte github repository (https://github.com/EarthByte/optAPM.git). Instructions on how to set up and run the TRM code, as well as the needed dependencies can be found in that repository. To run the TRM with the age error alterations as done in the SEIRF, I have added data folder that can be copied with the necessary files. It is therefore *not* necessary to copy this folder from the Earthbyte github.

- The SEIRF folder contains the files and Python codes needed to compute and plot the following:
    * The predicted motion path of Africa and true polar wander path with error
    * The predicted global trench kinematics and net lithosphere rotation with error
    * The predicted hotspot trails and hotspot source motion with error
    * The predicted reconstruction of kimberlite and large igneous province eruption sites

 The TRM with age errors produces a set of rotation files for each iteration and a csv file with the age range per iteration, which has to be copied into the optimisation folder in the SEIRF folder to compute the components of the SEIRF. The global trench kinematics, net lithosphere rotation, and reconstruction of kimberlite and large igneous province eruption sites take a long time to compute when done for each iteration, so their computation can be run on parallel cores or serially with the Calc_SEIRF_components file. This file is accompanied by the SEIRF_comp slurm file, which was created for the Utrecht University HPC. 

 The following Python files from the TRM code are necessary for the componenents of the SEIRF and have been copied to the SEIRF folder:
 - geoTools
 - net_rotation
 - optimised_rotation_updater
 - points_in_polygons
 - points_spatial_tree

The subduction_convergence_ptt file are a couple of functions copied from the PlateTectonicTools github repository (https://github.com/EarthByte/PlateTectonicTools.git) and adapted to compute the trench kinematics.

The SEIRF folder contains all rotation files for the Continent frame and Trench frame up to 1000 Ma. It is therefore not necessary to run the TRM with age error to compute parameters for these frames. The dependencies needed for the computing the SEIRF components and plotting the figures in the manuscript are the following:
- numpy
- pandas
- pmagpy
- mpi4py (if running in parallel)
- pygplates
- matplotlib
- cartopy
- future
