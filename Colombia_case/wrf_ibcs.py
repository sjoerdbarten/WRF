#wrf_ibcs.py
"""zero initial and boundary conditions in wrfinput and wrfbdy files for several tracers

Author: Ingrid Super, adapted from wrf_ct_tools.py from Denica Bozhinova, Michiel van der Molen and Wouter Peters
Last changes: 28-11-2014"""

"CO = 100 ppbv, O3 = 40 ppbv, NOx = 500 pptv"

import datetime as dtm
import os
import netCDF4 as nc
import numpy as np

#inpath = '/data/WRF/WRF3.2/Auke2/WRFV3/test/em_real'
inpath = '/home/WUR/barte035/WRFChem/WRFV3/run/'

zero_tracer_list = ['sulf',
                    'h2o2',
                    'ald',
                    'op1',
                    'op2',
                    'ora1',
                    'ora2',
                    'nh3',
                    'n2o5',
                    'no3',
                    'eth',
                    'ol2',
                    'olt',
                    'oli',
                    'tol',
                    'xyl',
                    'hono',
                    'hno4',
                    'ket',
                    'mgly',
                    'onit',
                    'csl',
                    'ho',
                    'ho2',
                    'hcl',
                    'ch3o2',
                    'ethp',
                    'ch3oh',
                    'c2h5oh',
                    'par',
                    'to2',
                    'cro',
                    'open',
                    'op3',
                    'c2o3',
                    'ro2',
                    'ano2',
                    'nap',
                    'xo2',
                    'xpar',
                    'isoprd',
                    'isopp',
                    'isopn',
                    'isopo2']

input_paths = [os.path.join(inpath,filename) for filename in os.listdir(inpath) if filename.startswith('wrfinput_')]
for each_path in input_paths:
    if os.path.exists(each_path):
        mf = nc.Dataset(each_path, mode='r+')
        bt = len(mf.dimensions['bottom_top'])
        sn = len(mf.dimensions['south_north'])
        we = len(mf.dimensions['west_east'])
        for each_tracer in zero_tracer_list:
            if each_tracer in mf.variables.keys():
               ncobj = mf.variables[each_tracer]
               tracer_ic = 1*[bt*[sn*[we*[0]]]]
               ncobj[:] = tracer_ic
            else:
                print 'Tracer %s skipped - not in the input file.'%each_tracer
                continue
        mf.close()
    else:
        print "Skipping %s - file path doesn't exist."%each_path
        continue

bdy_path = os.path.join(inpath,'wrfbdy_d01')
end_list = ['BXS', 'BXE', 'BYS', 'BYE', 'BTXS', 'BTXE', 'BTYS', 'BTYE']

if os.path.exists(bdy_path):
    mf = nc.Dataset(bdy_path, mode='r+')
    tm = len(mf.dimensions['Time'])
    bdyw = len(mf.dimensions['bdy_width'])
    bt = len(mf.dimensions['bottom_top'])
    for each_tracer in zero_tracer_list:
        for each_ending in end_list:
            new_tracer = '_'.join([each_tracer,each_ending])
            if new_tracer in mf.variables.keys():
                ncobj = mf.variables[new_tracer]
                len = ncobj.shape[3]
                tracer_bc = tm*[bdyw*[bt*[len*[0]]]]
                ncobj[:] = tracer_bc
            else:
                print 'Tracer %s skipped - not in the bdy file.'%new_tracer
                continue
    mf.close()
else:
    print "Skipping %s - file path doesn't exist."%bdy_path
    pass

