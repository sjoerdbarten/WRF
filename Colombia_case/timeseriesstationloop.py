import numpy as np
import matplotlib.pyplot as plt
import pylab
from scipy.stats import pearsonr

execfile("timeseriesstation.py")
WRFCO,OBSCO=timeseriesstation('CO','Bogota',0)
print(WRFCO.shape,OBSCO.shape)
print(WRFCO)
print(OBSCO)
print(OBSCO/WRFCO)

WRFO3,OBSO3=timeseriesstation('O3','Bogota',0)
WRFNO,OBSNO=timeseriesstation('NO','Bogota',0)

diffO3 = OBSO3-WRFO3
diffNO = OBSNO-WRFNO

print WRFO3.shape,OBSO3.shape,WRFNO.shape,OBSNO.shape
print pearsonr(OBSNO-WRFNO,OBSO3-WRFO3)

plt.plot(diffNO,'r-')
plt.plot(-diffO3,'b-')
plt.xlabel('Time')
plt.ylabel('Diff [ppb]')
plt.xlim([0, 24])
plt.ylim([0,70])
plt.show()
