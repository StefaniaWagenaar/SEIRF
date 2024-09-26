#--------------------
#   DECEMBER 2022
#   WRITTEN BY BRAM VAES
#
#   FUNCTIONS NEEDED FOR THE ROTATION OF POLES AND APWPS
#
#--------------------

import numpy as np 
import pandas as pd
import pmagpy.pmag as pmag
import pmagpy.pmagplotlib as pmagplotlib
import pmagpy.ipmag as ipmag
import pmagpy.contribution_builder as cb
from pmagpy import convert_2_magic as convert
import matplotlib.pyplot as plt

def sph2cart(lat,lon,r):
    """
    Converts an Euler pole into cartesian coordinates
    Euler pole should be given as (latitude, longitude, rotation angle)
    """
    x=r*np.cos(lat)*np.cos(lon)
    y=r*np.cos(lat)*np.sin(lon)
    z=r*np.sin(lat)
    return x, y, z

def rot2mat(pole):
    """
    Converts an Euler pole into a (3x3) rotation matrix
    Euler pole should be given as an array of length 3 (latitude, longitude, rotation angle)
    """ 
    gr=np.pi/180.
    E=sph2cart(pole[0]*gr,pole[1]*gr,1)   # converts pole to cartesian coordinates
    omega=pole[2]*gr                      # converts angle to radians
    R=np.zeros((3,3))                     # creates 3x3 matrix
    R[0][0]=E[0]*E[0]*(1-np.cos(omega)) + np.cos(omega)
    R[0][1]=E[0]*E[1]*(1-np.cos(omega)) - E[2]*np.sin(omega)
    R[0][2]=E[0]*E[2]*(1-np.cos(omega)) + E[1]*np.sin(omega)
    R[1][0]=E[1]*E[0]*(1-np.cos(omega)) + E[2]*np.sin(omega)
    R[1][1]=E[1]*E[1]*(1-np.cos(omega)) + np.cos(omega)
    R[1][2]=E[1]*E[2]*(1-np.cos(omega)) - E[0]*np.sin(omega)
    R[2][0]=E[2]*E[0]*(1-np.cos(omega)) - E[1]*np.sin(omega)
    R[2][1]=E[2]*E[1]*(1-np.cos(omega)) + E[0]*np.sin(omega)
    R[2][2]=E[2]*E[2]*(1-np.cos(omega)) + np.cos(omega)
    return R

def mat2rot(m):
    """
    Converts a (3x3) rotation matrix into an Euler pole
    Euler pole should be given as an array of length 3 (latitude, longitude, rotation angle)
    """
    gr=np.pi/180.
    lon = np.arctan2( m[0][2] - m[2][0] , m[2][1] - m[1][2] )
    term = np.sqrt( (m[2][1]-m[1][2])**2 + (m[0][2]-m[2][0])**2 + (m[1][0]-m[0][1])**2 )
    if term==0:
        pole=np.zeros(3)
    else:
        lat = np.arcsin( (m[1][0]-m[0][1]) / term )
        ang = np.arctan2( term , (m[0][0]+m[1][1]+m[2][2]-1) )
        pole = [lat/gr,lon/gr,ang/gr]
    return pole

def addpoles(p1,p2):
    """
    Adds two Euler poles
    Euler poles should be given as an array of length 3 (latitude, longitude, rotation angle)
    """
    A=rot2mat(p1)
    B=rot2mat(p2)
    C=np.dot(B,A)
    if np.array_equal(C,np.identity(3)) is True:
        result=np.zeros(3)
    else:
        result=mat2rot(C)
    return result

def stagepole(p_old,p_young):
    """
    Computes forward (total) stage pole
    For example, computes stage pole of Eurasia vs North America from 40 to 30 Ma
    Output is an Euler pole given as (latitude, longitude, rotation angle)
    Note that the rotation angle represents the total rotation during the stage
    """
    if np.array_equal(p_old,p_young) is True:
        stage_p=np.zeros(3)
    else:
        p_t2=p_old.copy()
        p_t2[2]=-1*p_t2[2]
        stage_p=addpoles(p_t2,p_young)
    return stage_p

def interpole(trp,sp,t2,t1,t):
    """
    Computes the total reconstruction pole at a desired age (t)
    Requires a total reconstruction pole at an older time (t2) and a stagepole describing the motion between t2 and t1
    Here, t1 is end of the motion stage described by the stage pole, such that t1 < t < t2
    Output is an Euler pole given as (latitude, longitude, rotation angle)
    """
    stagepole=sp.copy()
    delta=(t2-t)/(t2-t1)
    stagepole[2]=stagepole[2]*delta
    interp_pole=addpoles(trp,stagepole)
    return interp_pole 

def rot2frame(data,ID,t):
    """
    Computes the total reconstruction pole of a specific plate relative to the reference frame (ID=0)
    at a chosen age (t)

    Required input: a GPlates .rot file stored as a Numpy array 
    """
    plate_IDs=data[:,0]     # stores plate IDs in array
    ages=data[:,1]          # stores ages of rotation poles in array
    file_length=len(plate_IDs)  # computes length of rotation file
    
    trp=np.zeros(3)
    while (ID != 0):
        for j in range(file_length):              # loops through rotation file
            if (plate_IDs[j]==ID):
                if (ages[j]>=t):                  # finds age which is greater than t
                    pole_old=data[j,2:5]          # defines pole at t > t_interpolation
                    pole_young=data[j-1,2:5]      # defines pole at t < t_interpolation
                    pole_stage=stagepole(pole_old,pole_young)   # computes stage pole
                    pole_at_t=interpole(pole_old,pole_stage,ages[j],ages[j-1],t)    # computes pole at t_interpolation
                    trp=addpoles(trp,pole_at_t)   # computes new total pole
                    ID=data[j,5]                  # updates plate ID
    return trp


def pt_rot(EP, Lats, Lons):
    """
    Rotates points on a globe by an Euler pole rotation using method of
    Cox and Hart 1986, box 7-3.
    Parameters
    ----------
    EP : Euler pole list [lat,lon,angle] specifying the location of the pole;
    the angle is for a counterclockwise rotation about the pole
    Lats : list of latitudes of points to be rotated
    Lons : list of longitudes of points to be rotated
    Returns
    _________
    RLats : list of rotated latitudes
    RLons : list of rotated longitudes
    """
# gets user input of Rotation pole lat,long, omega for plate and converts
# to radians
    E = pmag.dir2cart([EP[1], EP[0], 1.])  # EP is pole lat,lon omega
    omega = np.radians(EP[2])  # convert to radians
    RLats, RLons = [], []
    for k in range(len(Lats)):
        if Lats[k] <= 90.:  # peel off delimiters
            # converts to rotation pole to cartesian coordinates
            A = pmag.dir2cart([Lons[k], Lats[k], 1.])
# defines cartesian coordinates of the pole A
            R = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]
            R[0][0] = E[0] * E[0] * (1 - np.cos(omega)) + np.cos(omega)
            R[0][1] = E[0] * E[1] * (1 - np.cos(omega)) - E[2] * np.sin(omega)
            R[0][2] = E[0] * E[2] * (1 - np.cos(omega)) + E[1] * np.sin(omega)
            R[1][0] = E[1] * E[0] * (1 - np.cos(omega)) + E[2] * np.sin(omega)
            R[1][1] = E[1] * E[1] * (1 - np.cos(omega)) + np.cos(omega)
            R[1][2] = E[1] * E[2] * (1 - np.cos(omega)) - E[0] * np.sin(omega)
            R[2][0] = E[2] * E[0] * (1 - np.cos(omega)) - E[1] * np.sin(omega)
            R[2][1] = E[2] * E[1] * (1 - np.cos(omega)) + E[0] * np.sin(omega)
            R[2][2] = E[2] * E[2] * (1 - np.cos(omega)) + np.cos(omega)
# sets up rotation matrix
            Ap = [0, 0, 0]
            for i in range(3):
                for j in range(3):
                    Ap[i] += R[i][j] * A[j]
# does the rotation
            Prot = pmag.cart2dir(Ap)
            RLats.append(Prot[1])
            RLons.append(Prot[0])
        else:  # preserve delimiters
            RLats.append(Lats[k])
            RLons.append(Lons[k])
    return RLons, RLats

def R_pole_space(ref_pole,ref_A95,obs_pole,obs_A95,local):
    """
    Compute rotation angle R with uncertainty dR using the pole-space approach of Ch.11 from Butler (1992)
    
    Required input:
    Observed pole -> [lon,lat]
    A95 of observed pole
    Reference pole -> [lon,lat]
    A95 of reference pole
    Sampling locality [lon,lat]

    Output:
    R, dR
    """

    p_r = pmag.angle(ref_pole,local) # compute angular distance p_r between reference pole and sampling locality
    p_o = pmag.angle(obs_pole,local) # compute angular distance p_o between observed pole and sampling locality
    s = pmag.angle(ref_pole,obs_pole) # compute angular distance s between observed pole and reference pole

    #print(p_r,p_o,s)

    # compute rotation angle R in degrees (eq. A72)
    R = np.rad2deg(np.arccos( (np.cos(np.deg2rad(s))-np.cos(np.deg2rad(p_o))*np.cos(np.deg2rad(p_r))) / (np.sin(np.deg2rad(p_o))*np.sin(np.deg2rad(p_r))) ))

    # compute dDx and dDo (eq. A74 and A75)
    delta_Dr = np.rad2deg (np.arcsin(np.sin(np.deg2rad(ref_A95))/np.sin(np.deg2rad(p_r))))
    delta_Do = np.rad2deg (np.arcsin(np.sin(np.deg2rad(obs_A95))/np.sin(np.deg2rad(p_o))))
    #print(delta_Do)

    # compute delta_R
    delta_R = 0.8*np.sqrt(delta_Dr**2+delta_Do**2)
    
    #print('R+dR = ',R[0],'+-',delta_R[0])

    return R[0], delta_R[0]