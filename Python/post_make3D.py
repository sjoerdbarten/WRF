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
lat_1 = 52.16 # 51.97
lon_1 = 5.74 #4.92

# for which date (initialization time) the wrfout needs to be analyzed.

(yy,mm,dd,HH) = (2018,6,5,00)

# which domain is selected

domain        = "d02"


# what is the name of the output file

outfilename = '3D.txt'
#--- User settings---#

infilename = '%s/wrfout_%s_%4d-%02d-%02d_%02d:00:00'%(datapath,domain,yy,mm,dd,HH)
print('opening %s ...'%infilename)

infile = nc.Dataset(infilename,'r')
print('Latitude_user = ',lat_1)
print('Longitude_user = ',lon_1)
xloc,yloc = wrf.to_np(wrf.ll_to_xy(infile,lat_1,lon_1))
print("xloc = "+str(xloc))
print("yloc = "+str(yloc))

print('Latitude_wrf:', infile.variables["XLAT"][0,yloc,xloc])
print('Longitude_wrf:', infile.variables["XLONG"][0,yloc,xloc])

nmodeltimesteps = infile.variables["U"][:,:,yloc,xloc].shape[0]
nmodellevels = infile.variables["U"][:,:,yloc,xloc].shape[1]

times = np.arange(0,nmodeltimesteps,1)
modellevels = np.arange(0,nmodellevels,1)

#We make an empty file which we later fill with the data
output = []
header = "IZ, Z(m),P(hPa),U(m/s),V(m/s),W(m/s),T(degC),Q(kg/kg),Qcloud(kg/kg)"
np.savetxt(outfilename,output,delimiter=',',header=header)
nodata = -999.000000

U = np.zeros((nmodeltimesteps,nmodellevels))
V = np.zeros((nmodeltimesteps,nmodellevels))
QVAPOR = np.zeros((nmodeltimesteps,nmodellevels))
QCLOUD = np.zeros((nmodeltimesteps,nmodellevels))
zx = np.zeros((nmodeltimesteps,nmodellevels))
pressx = np.zeros((nmodeltimesteps,nmodellevels))
Tx = np.zeros((nmodeltimesteps,nmodellevels))
Wx = np.zeros((nmodeltimesteps,nmodellevels))

for k in range(0,nmodeltimesteps,1):
        
        print('Processing timestep '+str(k)+' out of '+str(nmodeltimesteps-1))
	
        uvmet = wrf.getvar(infile,"uvmet",k).values
        
        U[k] = uvmet[0,:,yloc,xloc]
        V[k] = uvmet[1,:,yloc,xloc]
        QVAPOR[k] = infile.variables["QVAPOR"][k,:,yloc,xloc]
        QCLOUD[k] = infile.variables["QCLOUD"][k,:,yloc,xloc]
        zx[k] = wrf.getvar(infile,"z",k)[:,yloc,xloc] # grid point height
        pressx[k] = wrf.getvar(infile,"pressure",k)[:,yloc,xloc]
        Tx[k] = wrf.getvar(infile,"tc",k)[:,yloc,xloc]
        Wx[k] = wrf.getvar(infile,"wa",k)[:,yloc,xloc]	
	#Here we open the file and fill with the data
        outputfile = open(outfilename,"a")
        outputfile.write('***********************************************************************************************\n')
        output_emptyrow = np.vstack((times[k],nodata,nodata,nodata,nodata,nodata,nodata,nodata,nodata)).T
        np.savetxt(outputfile, output_emptyrow, fmt='%11.6f', delimiter=',', newline='\n')	
        output = np.vstack((modellevels,zx[k],pressx[k],U[k],V[k],Wx[k],Tx[k],QVAPOR[k],QCLOUD[k])).T
        np.savetxt(outputfile, output, fmt='%11.6f', delimiter=',', newline='\n')
        outputfile.close()
		
	
print('***Successfully finished make3D.py script***')
