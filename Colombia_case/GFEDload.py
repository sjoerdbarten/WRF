import h5py

GFED_data = '/archive/ESG/barte035/Colombia/GFED/fire_emissions_v4_R1_1293/data/GFED4.1s_2014.hdf5'

g = h5py.File(GFED_data, 'r')
xlat = g['lat'][:]
xlon = g['lon'][:]

DM_emissions = g['/emissions/01/DM'][:]

print(g)
print(xlat)
print(xlon)
print(DM_emissions)
