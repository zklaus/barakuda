#!/usr/bin/env python

#       B a r a K u d a
#
#       L. Brodeau, 2017]

import sys
import numpy as nmp
from PIL import Image
import string
import os
from netCDF4 import Dataset

narg = len(sys.argv)
if narg != 2:
    print 'Usage: '+sys.argv[0]+' <mesh_mask>'; sys.exit(0)
cf_mm = sys.argv[1]



cf_bmp = string.replace(os.path.basename(cf_mm), '.nc', '_orig.bmp')
cf_bmp = string.replace(os.path.basename(cf_bmp), '_orig.bmp4', '_orig.bmp')



# Opening mesh_mask:
f_mm = Dataset(cf_mm)
mask  = f_mm.variables['tmask'][0,0,:,:]
f_mm.close()


(nj, ni) = nmp.shape(mask)

print ' nj, ni =>', nj, ni

#imask= nmp.zeros((nj, ni), dtype=nmp.int8)

#imask[:,:] = mask[:,:]

#del mask

imask = (255*mask).astype(nmp.uint8)

# Then save it:
result = Image.fromarray(nmp.flipud(imask))
result.save(cf_bmp)
print ' *** Image '+cf_bmp+' saved!\n'
