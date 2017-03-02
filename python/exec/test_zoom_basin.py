#!/usr/bin/env python

#       B a r a K u d a
#
#       L. Brodeau, 2017]

import sys
import numpy as nmp
from netCDF4 import Dataset

narg = len(sys.argv)
if narg != 3:
    print 'Usage: '+sys.argv[0]+' <basin_mask> <basin_name>'; sys.exit(0)
cf_bm = sys.argv[1]
cname = sys.argv[2]





# Opening basin_mask:
f_bm = Dataset(cf_bm)
mask  = f_bm.variables['tmask'+cname][:,:]
f_bm.close()


(nj, ni) = nmp.shape(mask)

print ' nj, ni =>', nj, ni


(vjj , vji)  = nmp.where(mask==1)

#print vjj

nj_1 = max( nmp.min(vjj)-2 , 0    )
nj_2 = min( nmp.max(vjj)+2 , nj-1 )

ni_1 = max( nmp.min(vji)-2 , 0    )
ni_2 = min( nmp.max(vji)+2 , ni-1 )


#print " nmp.min(vji)-2 =", nmp.min(vjj)-2
#print " nmp.min(vjj)-2 =", nmp.min(vjj)-2


print ' nj_1,nj_2 =>', nj_1,nj_2
print ' ni_1,ni_2 =>', ni_1,ni_2
