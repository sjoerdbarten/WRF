import netCDF4 as nc
import scipy.interpolate
from pylab import *
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.axes as ax
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm

print 'This script makes a transect of vertical profiles over the whole output timeseries (or selected times). To use this use transectverticalprofile(wrffile,plotvar,levels,figsave) where wrffile is the path to the WRF-file, plotvar is any (4D) WRF-Chem variable, levels is how many levels (counting from surface) you want to plot. Figsave = 1 saves figure.'

def transectverticalprofile(wrffile,plotvar,levels,figsave):
	ncfile = nc.Dataset(wrffile,'r')
	
