#!/usr/bin/env python

# L. Brodeau, 2017

# Potential sigma0 density from potential temperature and salinity using
# function sigma0 (TEOS8) of barakuda_physics.py

import sys
#import os
import numpy as nmp
from netCDF4 import Dataset
from string  import replace

import barakuda_physics as bp

print '\n'

if len(sys.argv) != 4:
    print 'Usage: '+sys.argv[0]+' <NEMO grid_T file> <temperature_name> <salinity_name>'
    sys.exit(0)

cf_nemo_T = sys.argv[1]
cv_T      = sys.argv[2]
cv_S      = sys.argv[3]

cf_out = replace(cf_nemo_T, cf_nemo_T, 'sigma0_'+cf_nemo_T)

f_nemo_T = Dataset(cf_nemo_T)     # r+ => can read and write in the file... )
vtime   = f_nemo_T.variables['time_counter'][:] ; cu_time = f_nemo_T.variables['time_counter'].units
nav_lon = f_nemo_T.variables['nav_lon'][:,:]    ; cu_lon  = f_nemo_T.variables['nav_lon'].units
nav_lat = f_nemo_T.variables['nav_lat'][:,:]    ; cu_lat  = f_nemo_T.variables['nav_lat'].units
deptht  = f_nemo_T.variables['deptht'][:]       ; cu_dpt  = f_nemo_T.variables['deptht'].units
Nt = len(vtime)


for jt in range(Nt):

    print ' *** time record:', jt+1

    xtht  = f_nemo_T.variables[cv_T][jt,:,:,:]
    xsal  = f_nemo_T.variables[cv_S][jt,:,:,:]

    if jt == 0:
        (nk,nj,ni) = nmp.shape(xtht)
        xsg0 = nmp.zeros((nk,nj,ni))

    xsg0 = bp.sigma0(xtht, xsal) ; # computing Sigma0 at current time record !

    if jt == 0: 
        f_out = Dataset(cf_out, 'w', format='NETCDF3_CLASSIC')

        # Dimensions:
        f_out.createDimension('x', ni)
        f_out.createDimension('y', nj)
        f_out.createDimension('deptht', nk)
        f_out.createDimension('time_counter', None)

        id_lon  = f_out.createVariable('nav_lon','f4',('y','x',))  ; id_lon.units = cu_lon
        id_lat  = f_out.createVariable('nav_lat','f4',('y','x',))  ; id_lat.units = cu_lat
        id_dpt  = f_out.createVariable('deptht' ,'f4',('deptht',)) ; id_dpt.units = cu_dpt
        id_tim  = f_out.createVariable('time_counter' ,'f4',('time_counter',)) ; id_tim.units = cu_time

        id_lon[:,:] = nav_lon[:,:]
        id_lat[:,:] = nav_lat[:,:]
        id_dpt[:]   = deptht[:]

        id_sg0  = f_out.createVariable('sigma0','f4',('time_counter','deptht','y','x',))
        id_sg0.long_name = 'SIGMA0 density computed from T and S with TEOS8 / Jackett and McDougall (1994)'
        
        #f_out.About  = 'Bla bla'
        f_out.Author = 'barakuda [density_from_T_and_S.py] (https://github.com/brodeau/barakuda)'

    id_tim[jt]       = vtime[jt]
    id_sg0[jt,:,:,:] = xsg0[:,:,:]

f_out.close()

f_nemo_T.close()


print '\n'+cf_out+' sucessfully created!\n'
