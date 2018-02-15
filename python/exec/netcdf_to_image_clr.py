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

l_fake_coor = True
#l_fake_coor = False



narg = len(sys.argv)
if narg not in [4, 5]:
    print 'Usage: '+sys.argv[0]+' <netcdf_file.nc> <netcdf_variable> <image_extension (jpg,png,bmp,...)> (mutiple to field)'; sys.exit(0)

cf_nc = sys.argv[1]
cv_nc = sys.argv[2]
ciext = sys.argv[3]

imult = 1
if narg == 5: imult = int(sys.argv[4])


print imult
    
cfname, cncext = os.path.splitext(cf_nc)


cf_im = string.replace(os.path.basename(cf_nc), cncext, '.'+ciext)

print ' *** Will create image '+cf_im






# Reading data array:
f_nc = Dataset(cf_nc)
xfield  = imult*f_nc.variables[cv_nc][:,:]
f_nc.close()



(ny,nx) = nmp.shape(xfield)

#ifield = nmp.zeros((ny,nx), dtype=nmp.int16)

ifield = np.zeros((ny,nx,3), dtype=np.uint8)

print 'To be written!'
sys.exit(0)

xfield[:,:] = nmp.round(xfield[:,:], 0)

ia = xfield.astype(nmp.int16)

ifield[:,:,0] = ia*0.3


# Cleaning overshoots:
idx_too_small = nmp.where(ifield < 0)
ifield[idx_too_small] = 0
idx_too_large = nmp.where(ifield > 255)
ifield[idx_too_large] = 255

#print ifield[:,22]

ifield8 = ifield.astype(nmp.uint8)


image = Image.fromarray(nmp.flipud(ifield8))

# Then save it:
image.save(cf_im)
print ' *** Image '+cf_im+' saved!\n'
