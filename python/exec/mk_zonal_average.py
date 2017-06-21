#!/usr/bin/env python

# L. Brodeau, June 2012

import sys
from os.path import splitext,basename,dirname
import numpy as nmp
from netCDF4 import Dataset
from string import replace

import barakuda_tool as bt

# Defaults:
cv_lon = 'lon'
cv_lon = 'lat'

narg = len(sys.argv)
if not narg in [ 5 , 7]:
    print 'Usage: '+sys.argv[0]+' <FILE_lat-lon.nc> <variable> <mesh_mask.nc> <mask_name> (<name_longitude> <name_latitude>)'
    sys.exit(0)

cf_in  = sys.argv[1]
cv_in  = sys.argv[2]

cf_msk  = sys.argv[3]
cv_msk  = sys.argv[4]

if narg == 7:
    cv_lon = sys.argv[5]
    cv_lat = sys.argv[6]

print ''

cn_file, cn_ext = splitext(cf_in)
cn_file = replace(cn_file, cv_in+'_', '')
cpath='.'
if dirname(cf_in) != '': cpath=dirname(cf_in)
cf_out = cpath+'/zonal_'+cv_in+'_'+basename(cn_file)+'.nc'
#print ' *[mk_zonal_average.py]* file to write: ', cf_out ; sys.exit(0)

bt.chck4f(cf_msk)
f_msk = Dataset(cf_msk)
Ndim = len(f_msk.variables[cv_msk].dimensions)
if   Ndim == 4:
    xmsk = f_msk.variables[cv_msk][0,0,:,:]
elif Ndim == 3:
    xmsk = f_msk.variables[cv_msk][0,:,:]
elif Ndim == 2:
    xmsk = f_msk.variables[cv_msk][:,:]
else:
    print ' ERROR (mk_zonal_average.py) => weird shape for your mask array!'
    sys.exit(0)    
f_msk.close()

bt.chck4f(cf_in)
f_in = Dataset(cf_in)

list_var = f_in.variables.keys()

Ndim = len(f_in.variables[cv_lon].dimensions)

cunt_lon = f_in.variables[cv_lon].units
cunt_lat = f_in.variables[cv_lat].units

if Ndim == 1:
    # Extracting the longitude and 1D array:
    vlon     = f_in.variables[cv_lon][:]
    # Extracting the latitude 1D array:
    vlat     = f_in.variables[cv_lat][:]

elif Ndim == 2:
    # We suppose it is NEMO nav_lon and nav_lat...
    xlon = f_in.variables[cv_lon][:,:]
    xlat = f_in.variables[cv_lat][:,:]
    (nj0,ni0) = nmp.shape(xlon)
    vlon = nmp.zeros(ni0)
    vlat = nmp.zeros(nj0)
    vlon[:] = xlon[nj0/8,:]
    ji_lat0 = nmp.argmax(xlat[nj0-1,:])
    vlat[:] = xlat[:,ji_lat0]
    del xlon, xlat, nj0, ni0
else:
    print ' ERROR (mk_zonal_average.py) => weird shape for your longitude array!'
    sys.exit(0)

for cvt in [ 'time', 'time_counter' ]: 
    if cvt in list_var:
        vtime     = f_in.variables[cvt][:]
        cunt_time = f_in.variables[cvt].units
#print 'TIME: ', cunt_time, '\n'

# Field !!!
xfield      = f_in.variables[cv_in][:,:,:]
cunit_field = f_in.variables[cv_in].units
clgnm_field = f_in.variables[cv_in].long_name

f_in.close()


Nt = len(vtime)
#print ' Nt = '+str(Nt)


# Checking dimensions
# ~~~~~~~~~~~~~~~~~~~
[ Nt, nj, ni ] = xfield.shape
print ' *[mk_zonal_average.py]* dimension of "'+cv_in+'" => ', ni, nj, Nt

Fzonal = bt.mk_zonal(xfield, xmsk)

# Output netCDF file:
#######################
f_out = Dataset(cf_out, 'w',format='NETCDF3_CLASSIC')

# Dimensions:
f_out.createDimension('lat', nj)
f_out.createDimension('time', None)

# Variables
id_lat = f_out.createVariable('lat','f4',('lat',))
id_tim = f_out.createVariable('time','f4',('time',))

id_f1  = f_out.createVariable(cv_in,'f4',('time','lat',))
id_f1.units     = cunit_field
id_f1.long_name = clgnm_field

id_f2  = f_out.createVariable(cv_in+'_mean','f4',('lat',))
id_f2.units     = cunit_field
id_f2.long_name = clgnm_field

id_f3  = f_out.createVariable(cv_in+'_anom','f4',('time','lat',))

f_out.about = 'Diagnostics created with BaraKuda (https://github.com/brodeau/barakuda)'

# Filling variables:
id_lat[:] = vlat[:]

Z_time_mean = nmp.mean(Fzonal[:,:], axis=0)

for jt in range(Nt):
    id_tim[jt] = vtime[jt]
    id_f1[jt,:] = Fzonal[jt,:]
    id_f3[jt,:] = Fzonal[jt,:] - Z_time_mean

id_f2[:] = Z_time_mean[:]

f_out.close()
        

print ' *[mk_zonal_average.py]* Wrote file '+cf_out+' !\n'

