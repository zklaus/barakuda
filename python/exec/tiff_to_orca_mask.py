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
import datetime

narg = len(sys.argv)
if narg != 2:
    print 'Usage: '+sys.argv[0]+' <mesh_mask>'; sys.exit(0)
cf_mm = sys.argv[1]

cim_pac = 'mask_pac.tiff'
cim_atl = 'mask_atl.tiff'
cim_ind = 'mask_ind.tiff'

for cf in [ cim_pac , cim_atl , cim_ind ]:
    if not os.path.exists(cf):
        print 'PROBLEM! image '+cf+' is not here!!!' ; sys.exit(0)
    else:
        print ' *** good, '+cf+' is here...'


# Opening mesh_mask:
f_mm = Dataset(cf_mm)
nav_lon = f_mm.variables['nav_lon'][:,:]
nav_lat = f_mm.variables['nav_lat'][:,:]
mask    = f_mm.variables['tmask'][0,0,:,:]
f_mm.close()

(nj,ni) = nmp.shape(nav_lon)


cf_bm = string.replace(os.path.basename(cf_mm), 'mesh_', 'basin_')


# Opening Images:
pic_pac = Image.open(cim_pac)
pic_atl = Image.open(cim_atl)
pic_ind = Image.open(cim_ind)

im_pac_array = nmp.flipud( nmp.array(pic_pac) )
im_atl_array = nmp.flipud( nmp.array(pic_atl) )
im_ind_array = nmp.flipud( nmp.array(pic_ind) )

(njp,nip) = im_pac_array.shape
(nja,nia) = im_atl_array.shape
(nji,nii) = im_ind_array.shape

if (njp,nip) != (nj,ni) or (nja,nia) != (nj,ni) or (nji,nii) != (nj,ni):
    print 'ERRO: something is wrong with the shapes of 2D arrays:'
    print nj,ni
    print njp,nip
    print nja,nia
    print nji,nii
    sys.exit(0)


mask_pac = nmp.zeros((nj,ni))
idx_pac = nmp.where(im_pac_array > 0)
mask_pac[idx_pac] = 1

mask_atl = nmp.zeros((nj,ni))
idx_atl = nmp.where(im_atl_array > 0)
mask_atl[idx_atl] = 1

mask_ind = nmp.zeros((nj,ni))
idx_ind = nmp.where(im_ind_array > 0)
mask_ind[idx_ind] = 1

now = datetime.datetime.now()
cdate = now.strftime("%Y-%m-%d")



f_out = Dataset(cf_bm, 'w', format='NETCDF3_CLASSIC')

# Dimensions:
f_out.createDimension('x', ni)
f_out.createDimension('y', nj)


id_lon = f_out.createVariable('nav_lon' ,'f4',('y','x',))
id_lat = f_out.createVariable('nav_lat' ,'f4',('y','x',))

id_pac  = f_out.createVariable('tmaskpac','f4',('y','x',))
id_atl  = f_out.createVariable('tmaskatl','f4',('y','x',))
id_ind  = f_out.createVariable('tmaskind','f4',('y','x',))


id_lon[:,:] = nav_lon[:,:]
id_lat[:,:] = nav_lat[:,:]

id_pac[:,:] = mask_pac[:,:]*mask[:,:]
id_atl[:,:] = mask_atl[:,:]*mask[:,:]
id_ind[:,:] = mask_ind[:,:]*mask[:,:]

f_out.About  = 'ORCA025, mask for main ocean basins, created with orca_mesh_mask_to_bitmap.py, Gimp, and tiff_to_orca_mask.py, '+cdate+'.'
f_out.Author = 'L. Brodeau (https://github.com/brodeau/barakuda)'

f_out.close()



print cf_bm+' created!!!'

