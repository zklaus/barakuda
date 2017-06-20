#!/usr/bin/env python

# L. Brodeau, June 2012

import sys
from os.path import splitext,basename,dirname
import numpy as nmp
from netCDF4 import Dataset
from string import replace

import barakuda_tool as bt

rmv    = -9999.

# Defaults:
cv_lon = 'lon'
cv_lon = 'lat'

narg = len(sys.argv)
if not narg in [ 3 , 5]:
    print 'Usage: '+sys.argv[0]+' <FILE_lat-lon.nc> <variable> (<name_longitude> <name_latitude>)'
    sys.exit(0)

cf_in  = sys.argv[1]
cv_in  = sys.argv[2]

if narg == 5:
    cv_lon = sys.argv[3]
    cv_lat = sys.argv[4]

cn_file, cn_ext = splitext(cf_in)
cn_file = replace(cn_file, cv_in+'_', '')
cf_out = dirname(cf_in)+'/zonal_'+cv_in+'_'+basename(cn_file)+'.nc'



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

print 'LONGITUDE: ', cunt_lon
print 'LATITUDE: ', cunt_lat



for cvt in [ 'time', 'time_counter' ]: 
    if cvt in list_var:
        vtime     = f_in.variables[cvt][:]
        cunt_time = f_in.variables[cvt].units
print 'TIME: ', cunt_time, '\n'

# Field !!!
#rmv    = f_in.variables[cv_in]._FillValue
xfield      = f_in.variables[cv_in][:,:,:]
cunit_field = f_in.variables[cv_in].units
clgnm_field = f_in.variables[cv_in].long_name

#print 'Missing value for '+cv_in+' is : ', rmv, '\n'

f_in.close()


Nt = len(vtime)
print ' Nt = '+str(Nt)




# Checking dimensions
# ~~~~~~~~~~~~~~~~~~~
[ Nt, nj, ni ] = xfield.shape
print ' DIMENSION =>  ni, nj, Nt = ', ni, nj, Nt


VZ = nmp.zeros((Nt,nj))

for jt in range(Nt):

    cjt  = '%3.3d' %(jt+1)

    # Zonally-averaging:
    for jj in range(nj):
    
        cpt = 0
    
        for ji in range(ni):
            val = xfield[jt,jj,ji]
            if val != rmv:
                cpt = cpt + 1
                VZ[jt,jj] = VZ[jt,jj] + val

        if cpt >= 1 : VZ[jt,jj] = VZ[jt,jj]/cpt



print '\n Creating netcdf file! ()'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



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

Z_time_mean = nmp.mean(VZ[:,:], axis=0)

for jt in range(Nt):
    id_tim[jt] = vtime[jt]
    id_f1[jt,:] = VZ[jt,:]
    id_f3[jt,:] = VZ[jt,:] - Z_time_mean

id_f2[:] = Z_time_mean[:]

f_out.close()
        

print '\n *** Wrote file '+cf_out+' !\n'

