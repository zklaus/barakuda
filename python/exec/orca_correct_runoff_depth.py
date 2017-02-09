#!/usr/bin/env python


# L. Brodeau, 2017

import sys
import os
import numpy as nmp
from netCDF4 import Dataset
from string import replace

cv_rnf_dept = 'rodepth'
cv_bathy_m  = 'Bathymetry'

if len(sys.argv) != 4:
    print 'Usage: '+sys.argv[0]+' <runoff_depth.nc> <bathy_meter.nc> <min depth (m)>'
    sys.exit(0)
cf_rnf_dept  = sys.argv[1]
cf_bathy_m   = sys.argv[2]
cmin_dept    = sys.argv[3] ; rmin_depth = float(cmin_dept)

cf_new = replace(cf_rnf_dept, '.nc', '_min'+cmin_dept+'.nc')

os.system('rm -f '+cf_new)
os.system('cp '+cf_rnf_dept+' '+cf_new)

print '\n'


f_bathy = Dataset(cf_bathy_m)
xbathy = f_bathy.variables[cv_bathy_m][:,:]
f_bathy.close()

(Nj,Ni) = nmp.shape(xbathy)

xnew = nmp.zeros((Nj,Ni))

print '\n'


# Opening the Netcdf file:
f_new = Dataset(cf_new, 'r+')     # r+ => can read and write in the file... )
print 'File ', cf_new, 'is open...\n'

# Extracting tmask at surface level:
xtemp  = f_new.variables[cv_rnf_dept][:,:,:]

xnew[:,:] = xtemp[0,:,:]
idx = nmp.where( (xbathy[:,:] >= rmin_depth) & (xnew[:,:] < rmin_depth) )

#print idx

xnew[idx] = rmin_depth

# Updating:
f_new.variables[cv_rnf_dept][0,:,:] = xnew[:,:]

f_new.Author = 'L. Brodeau (orca_correct_runoff_depth.py of Barakuda)'

f_new.close()

print cf_new+' sucessfully created!'

