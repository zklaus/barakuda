#!/usr/bin/env python

#       B a r a K u d a
#
#       L. Brodeau, 2017]

import sys
import numpy as nmp
import string
import os
from netCDF4 import Dataset

import barakuda_orca as bo


narg = len(sys.argv)
if narg != 4:
    print 'Usage: '+sys.argv[0]+' <mesh_mask> <lon> <lat>'; sys.exit(0)
cf_mm = sys.argv[1]
clon = sys.argv[2] ; rlon = float(clon)
clat = sys.argv[3] ; rlat = float(clat)

print rlon, rlat


# Opening mesh_mask:
f_mm = Dataset(cf_mm)
nav_lon = f_mm.variables['nav_lon'][:,:]
nav_lat = f_mm.variables['nav_lat'][:,:]
#mask    = f_mm.variables['tmask'][0,0,:,:]
f_mm.close()

(nj,ni) = nmp.shape(nav_lon)



(ji, jj) = bo.ij_from_xy(rlon, rlat, nav_lon, nav_lat)


print 'RESULT => ji, jj =', ji,jj
print '  lon, lat =', nav_lon[jj,ji], nav_lat[jj,ji], '\n'

