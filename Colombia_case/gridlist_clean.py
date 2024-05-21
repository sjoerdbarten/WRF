#!/usr/local/bin/python
from pylab import *
from numpy import *
#from pyhdf import SD
#from Scientific.IO.NetCDF import *
from netCDF4 import *
import sys
import subprocess
import os
from scipy import stats
import itertools
from collections import OrderedDict
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.ticker as mtick
#import h5py
from mpl_toolkits.basemap import Basemap

##########################################################################################
#                        Program to re-grid satellite data                               #
#                                                                                        #
#                            October  2017                                               #  
#                      A. Lorente based on gridlist by KFB                               #
##########################################################################################
#
# Aim of this code is to re-grid satellite data into a 0.25X0.25 grid by averaging
#
##########################################################################################




def gridlist(lon4, lat4, im, jm):#, nind, iind, jind):
# Function gridlist:
#    !--------------------------------------------------------------  
#    ! 
#    ! Given the integer corner positions (ilat,lon), 
#    ! return a list of grid cell indices (iind,jind), dimension nind
#    !
#    ! Change - 23/05/2011 Folkert Boersma and Joost Maasakkers
#    !    Reverse order of maxm and maxk for OMI, and define maxm
#    !    dynamically.
#    !
#    ! Routine requires:
#    !   - That (1,4) and (2,3) are opposite edges
#    !   - That maxm > maxk
#    !   - That { pair (1,3) have highest latitude, (2,4) lowest } or
#    !          { pair (2.4) have highest latitude, (1,3) lowest }
#    !
#    !     1 - - - - - - - - 3
#    !     |                 |
#    !     |                 |
#    !     2 - - - - - - - - 4     (maxk = 2, maxm = 8)
#    !
#    !               Folkert Boersma, Harvard University, March 2006
#    !                                    Henk Eskes, KNMI, Aug 2017
#    !                                    Alba Lorente, WUR, Oct 2017
#    !                 Translate to python from original fortran code
#    !                 and adaptation to GOME-2 data                                           
#    !--------------------------------------------------------------


#lon4,lat4: contains the corner coordinates of each pixel

#!Number of sampling points along the obs. ground pixel
    #!integer, parameter     ::  maxk = 2   ! settings Folkert
    #!integer, parameter     ::  maxm = 16
    #!integer, parameter     ::  maxk = 4
    #!integer, parameter     ::  maxm = 32
    maxk = 4   
    maxm = 32

    iind=[];jind=[]
    xcorner = np.empty([4])
    ycorner = np.empty([4])
    xcorner[0]=lon4[3]; ycorner[0]=lat4[3]
    xcorner[1]=lon4[2]; ycorner[1]=lat4[2]
    xcorner[2]=lon4[1]; ycorner[2]=lat4[1]
    xcorner[3]=lon4[0]; ycorner[3]=lat4[0]

# Longitudes from 0-360:
    for i4 in range(4):
        if (xcorner[i4] < 0): xcorner[i4] = xcorner[i4] + 360.0
    
    d1 = abs(xcorner[0]-xcorner[1])
    d2 = abs(xcorner[0]-xcorner[2])
    d3 = abs(xcorner[0]-xcorner[3])
    if ((d1 > 180.0) or (d2 > 180.0) or (d3 > 180.0)):
        for i4 in range(4):
            if (xcorner[i4] >180.): xcorner[i4] = xcorner[i4] - 360.0
                
    # sub-cell method
    i = 0
    for k in range(1,maxk+1,1):
        for m in range(1,maxm+1,1):
            # longitude of sub-cell centre of groundpixel	
            # scan simultaneously from 1 -> 3 and 2 -> 4
            lon_low  = xcorner[0] + \
               (xcorner[2]-xcorner[0])*(2*m-1)/(2.*maxm)
            lon_high = xcorner[1] + \
               (xcorner[3]-xcorner[1])*(2*m-1)/(2.*maxm)
            lon_cell = lon_low + (lon_high-lon_low)*(2*k-1)/(2.*maxk)
            # latitude of sub-cell centre of groundpixel	
            lat_low  = ycorner[0] + \
               (ycorner[2]-ycorner[0])*(2*m-1)/(2.*maxm)
            lat_high = ycorner[1] + \
               (ycorner[3]-ycorner[1])*(2*m-1)/(2.*maxm)
            lat_cell = lat_low+(lat_high-lat_low)*(2*k-1)/(2.*maxk)
            
            # calculating the indices of the  grid-cell
            x = 0.5 + ( lon_cell + 180.0 ) * im / 360.0
            y = 0.5 + ( lat_cell +  90.0 ) * jm / 180.0
            i = i+1
            xi = ( (round(x) - 1) % im ) + 1
            yi = ( (round(y) - 1) % jm ) + 1
            iind.append(int(xi))
            jind.append(int(yi))

    nind = i
    iind = np.array(iind);jind = np.array(jind)
    return nind,iind,jind
