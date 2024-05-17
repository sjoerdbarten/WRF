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

outfilename = '2D.txt'

#--- User settings---#

#--- This function calculates the boundary layer height and is consistent with all boundary layer schemes ---#
def pbldepth(u,v,z,T,nz):
	ricrit=0.2
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
print('Latitude_user = ',lat_1)
print('Longitude_user = ',lon_1)
xloc,yloc = wrf.to_np(wrf.ll_to_xy(infile,lat_1,lon_1))
print("xloc = "+str(xloc))
print("yloc = "+str(yloc))
print('Latitude_wrf:', wrf.getvar(infile,"lat",0)[yloc,xloc].values)
print('Longitude_wrf:',wrf.getvar(infile,"lon",0)[yloc,xloc].values)
nmodeltimesteps = infile.variables["U"][:,:,yloc,xloc].shape[0]
nmodellevels = infile.variables["U"][:,:,yloc,xloc].shape[1]

times = np.arange(0,nmodeltimesteps,1)

uvmet = wrf.getvar(infile,"uvmet",timeidx=wrf.ALL_TIMES)
U = uvmet.values[0,:,:,yloc,xloc]
V = uvmet.values[1,:,:,yloc,xloc]
ht = wrf.getvar(infile,"height",wrf.ALL_TIMES)
z1 = ht.values[:,:,yloc,xloc]

pbl2 = np.zeros(nmodeltimesteps)
for k in range(0,nmodeltimesteps,1):
	th = wrf.getvar(infile,"theta",timeidx=k)
	T1P = th.values[:,yloc,xloc]
	pbl2[k] = pbldepth(U[k],V[k],z1[k],T1P,nmodellevels)

pbl1 = infile.variables["PBLH"][:,yloc,xloc]
HFX = infile.variables["HFX"][:,yloc,xloc]
TSK = infile.variables["TSK"][:,yloc,xloc]-273.15
GRDFLX = infile.variables["GRDFLX"][:,yloc,xloc]
LH = infile.variables["LH"][:,yloc,xloc]
UST = infile.variables["UST"][:,yloc,xloc]
T2 = infile.variables["T2"][:,yloc,xloc]-273.15

uvmet10 = wrf.getvar(infile,"uvmet10",timeidx=wrf.ALL_TIMES)
U10 = uvmet10.values[0,:,yloc,xloc]
V10 = uvmet10.values[1,:,yloc,xloc]

SWDOWN = infile.variables["SWDOWN"][:,yloc,xloc]
ALBEDO = infile.variables["ALBEDO"][:,yloc,xloc]
GLW = infile.variables["GLW"][:,yloc,xloc]
GSW = SWDOWN*(1.-ALBEDO)
LU_INDEX = infile.variables["LU_INDEX"][:,yloc,xloc]
Q2 = infile.variables["Q2"][:,yloc,xloc]
PSFC = infile.variables["PSFC"][:,yloc,xloc]/100.
SMOIS = infile.variables["SMOIS"][:,0,yloc,xloc]
TSLB = infile.variables["TSLB"][:,0,yloc,xloc]
EMISS = infile.variables["EMISS"][:,yloc,xloc]
XLAT = infile.variables["XLAT"][:,yloc,xloc]
XLONG = infile.variables["XLONG"][:,yloc,xloc]

output = np.vstack((times,PSFC,T2,Q2,U10,V10,UST,SWDOWN,GSW,GLW,TSK,HFX,LH,GRDFLX,pbl1,pbl2,LU_INDEX,TSLB,SMOIS,ALBEDO,EMISS,XLAT,XLONG)).T
header = "Time,P_surf(hPa),T2m(oC),Q2m(kg/kg),U10(m/s),V10(m/s),Ustar(m/s),SWdown(W/m2),Net_SW(W/m2),GLW(W/m2),Tskin(oC),H(W/m2),LvE(W/m2),G(W/m2),pbl_1(m),pbl_2(m),luindex(-),soiltemp(K),soil_moisture(-),albedo(-),emissivity(-),longitude,latitude"

ftmx= "%12.4f"

np.savetxt(outfilename,output,delimiter=',',header=header,fmt=ftmx)

print('***Successfully finished make2D.py script***')
