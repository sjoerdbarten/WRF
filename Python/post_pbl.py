import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import wrf
np.bool = np.bool_

#--- User settings---#

home = os.path.expanduser('~') # set home directory
datapath = '%s/atmmod/WRF/run'%home # set full datapath

# for which point (latitude,longitude) the pbl is computed
lat_1 = 51.97
lon_1 = 4.92

# for which date (initialization time) the wrfout needs to be analyzed.

(yy,mm,dd,HH) = (2018,6,5,00)

# which domain is selected

domain        = "d02"

# what is the name of the output file

outfilename = 'pblheight.txt'

#--- User settings---#

#--- This function calculates the boundary layer height and is consistent with all boundary layer schemes ---#
def pbldepth(u,v,z,T,nz):
	ricrit=0.25
	ribulk = np.zeros(nz)
	for i in range(0,nz,1):
		ribulk[i] = 9.81/300.*(z[i]-z[0])*(T[i]-T[0])/(u[i]**2+v[i]**2)
	i = 0
	while ribulk[i] < ricrit:
		i = i+1
	pblhheight = z[i-1]+(z[i]-z[i-1])/(ribulk[i]-ribulk[i-1])*(ricrit-ribulk[i-1])
	return pblhheight
#--- This function calculates the boundary layer height and is consistent with all boundary layer schemes ---#

infilename = '%s/wrfout_%s_%4d-%02d-%02d_%02d:00:00'%(datapath,domain,yy,mm,dd,HH)
print('opening %s ...'%infilename)

infile = nc.Dataset(infilename,'r')

xloc,yloc = wrf.to_np(wrf.ll_to_xy(infile,lat_1,lon_1))

print('Latitude:', infile.variables["XLAT"][0,yloc,xloc])
print('Longitude:', infile.variables["XLONG"][0,yloc,xloc])

nmodeltimesteps = infile.variables["U"][:,:,yloc,xloc].shape[0]
nmodellevels = infile.variables["U"][:,:,yloc,xloc].shape[1]

times = np.arange(0,nmodeltimesteps,1)
pblh = np.zeros(nmodeltimesteps)
for k in range(0,nmodeltimesteps,1):
	p1 = wrf.getvar(infile,"pres",k)[:,yloc,xloc]
	z = wrf.getvar(infile,"z",k)[:,yloc,xloc]
	qv1 = wrf.getvar(infile,"QVAPOR",k)[:,yloc,xloc]*1000.
	tc1 = wrf.getvar(infile,"tc",k)[:,yloc,xloc]
	u1 = wrf.getvar(infile,"ua",k)[:,yloc,xloc]
	v1 = wrf.getvar(infile,"va",k)[:,yloc,xloc]
	theta1 = wrf.getvar(infile,"theta",timeidx=k).values[:,yloc,xloc]
	pblh[k] = pbldepth(u1,v1,z,theta1,nmodellevels)
	

	
output = np.vstack((times,pblh)).T
header = "Time,PBLheight(m), lat: "+str(lat_1)+",lon: " +str(lon_1)
np.savetxt(outfilename,output,delimiter=',',header=header,fmt='%12.4f')

print('***Successfully finished pbl.py script***')

