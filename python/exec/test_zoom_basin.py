#!/usr/bin/env python

#       B a r a K u d a
#
#       L. Brodeau, 2017]

import sys
import numpy as nmp
from netCDF4 import Dataset

import barakuda_orca as bo

narg = len(sys.argv)
if narg != 3:
    print 'Usage: '+sys.argv[0]+' <basin_mask> <basin_name>'; sys.exit(0)
cf_bm = sys.argv[1]
cname = sys.argv[2]





# Opening basin_mask:
f_bm = Dataset(cf_bm)
mask  = f_bm.variables['tmask'+cname][:,:]
f_bm.close()


(i1,j1,i2,j2) = bo.shrink_domain(mask)


print ' i1, i2 =>', i1, i2
print ' j1, j2 =>', j1, j2

