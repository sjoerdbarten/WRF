import netCDF4 as nc
import csv
import pandas as pd
import numpy as np

wrffile = '/archive/ESG/barte035/WRFV3/wrfout_d01_2014-01-01_00:00:00'
data = nc.Dataset(wrffile,'r')

varlist=['PSFC','P','U','V','W','QVAPOR','T','o3','no','no2','co']
df = pd.DataFrame(index=np.arange(0,data.variables['U'].shape[0],1),columns=np.arange(0,((len(varlist)-1)*data.variables['U'].shape[1]+1),1))
print(df)

for v in range(0,len(varlist)):
	for z in range(0,data.variables['U'].shape[1]):
		for t in range(0,data.variables['U'].shape[0]):
			if varlist[v] in ['PSFC']:
				df.iloc[t,0] = data.variables['PSFC'][t,48,34]
			if varlist[v] in ['P']:
				df.iloc[t,z+((v-1)*data.variables['U'].shape[1])+1] = data.variables['P'][t,z,48,34]+data.variables['PB'][t,z,48,34]
			if varlist[v] in ['U','V','W','QVAPOR','o3','no','no2','co']:
				df.iloc[t,z+((v-1)*data.variables['U'].shape[1])+1] = data.variables[varlist[v]][t,z,48,34]
			if varlist[v] in ['T']:
				df.iloc[t,z+((v-1)*data.variables['U'].shape[1])+1] = (data.variables['T'][t,z,48,34] + 300.) * (((data.variables['P'][t,z,48,34]+data.variables['PB'][t,z,48,34]) / 100000.) ** 0.2854)
			if varlist[v] not in ['PSFC']:
				df.rename(columns={z+((v-1)*data.variables['U'].shape[1])+1:str(varlist[v])+str(z)},inplace=True)
			if varlist[v] in ['PSFC']:
				df.rename(columns={0:'PSFC'},inplace=True)			
		print(z)
	print(df)
	
df.to_csv('/archive/ESG/barte035/Colombia/1DinputBogota_includingchem.txt')

print(df.columns)
print(df[['P0','P1','P2']])
print(df[['U0','U1','U2']])
print(df[['V0','V1','V2']])
print(df[['W0','W1','W2']])
print(df[['QVAPOR0','QVAPOR1','QVAPOR2']])
print(df[['T0','T1','T2']])
print(df['PSFC'])
print(df[['o30','o31','o32']])
print(df[['no0','no1','no2']])
print(df[['no20','no21','no22']])
print(df[['co0','co1','co2']])

#u = data.variables['U'][:,:,48,34] #m s-1
#v = data.variables['V'][:,:,48,34] #m s-1
#q = data.variables['QVAPOR'][:,:,48,34]
#p = data.variables['P'][:,:,48,34]+data.variables['PB'][:,:,48,34]
#T = (data.variables['T'][:,:,48,34] + 300.) * ((p[:,:] / 100000.) ** 0.2854)


#print(u)
#print(v)
#print(q)
#print(p)
#print(T)
