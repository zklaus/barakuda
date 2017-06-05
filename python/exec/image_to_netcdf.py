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

l_fake_coor = True
#l_fake_coor = False



narg = len(sys.argv)
if narg != 2:
    print 'Usage: '+sys.argv[0]+' <image>'; sys.exit(0)

cf_im = sys.argv[1]

cfname, cfext = os.path.splitext(cf_im)


#(nj,ni) = nmp.shape(nav_lon)

cf_nc = string.replace(os.path.basename(cf_im), cfext, '.nc')

# Opening Images:
print ' *** Opening image '+cf_nc

pic = Image.open(cf_im)


print nmp.shape(pic)

#sys.exit(0)
#xtmp = nmp.flipud( nmp.array(pic) )

(ny,nx,nrgb) = nmp.shape(pic)

xpic = nmp.array(pic)

if nrgb != 3: print ' ERROR: we expect 3 colors! RGB!'; sys.exit(0)


if l_fake_coor:
    # Prepare coordinates if needed:
    vlon = nmp.zeros(nx) ; dx = 360./float(nx)
    for ji in range(nx): vlon[ji] = (float(ji) + 0.5)*dx
    
    vlat = nmp.zeros(ny) ; dy = 180./float(ny)
    for jj in range(ny): vlat[jj] = -90 + (float(jj) + 0.5)*dy
    #print vlat[:]
    #sys.exit(0)


f_out = Dataset(cf_nc, 'w', format='NETCDF4')

# Dimensions:

cdim_x = 'x'
cdim_y = 'y'
if l_fake_coor:
    cdim_x = 'lon'
    cdim_y = 'lat'
f_out.createDimension(cdim_x, nx)
f_out.createDimension(cdim_y, ny)

if l_fake_coor:
    id_lon  = f_out.createVariable('lon' ,'f4',(cdim_x,))
    id_lat  = f_out.createVariable('lat' ,'f4',(cdim_y,))
    id_lon[:] = vlon[:]
    id_lat[:] = vlat[:]




id_red  = f_out.createVariable('red','f4',(cdim_y,cdim_x,))
id_red.long_name = 'Red (of RGB)'

id_green  = f_out.createVariable('green','f4',(cdim_y,cdim_x,))
id_green.long_name = 'Green (of RGB)'

id_blue  = f_out.createVariable('blue','f4',(cdim_y,cdim_x,))
id_blue.long_name = 'Blue (of RGB)'

id_red[:,:]   = nmp.flipud(xpic[:,:,0])
id_green[:,:] = nmp.flipud(xpic[:,:,1])
id_blue[:,:]  = nmp.flipud(xpic[:,:,2])

f_out.About  = 'Image '+cf_im+' converted to netcdf.'
f_out.Author = 'Generated with image_to_netcdf.py of BARAKUDA (https://github.com/brodeau/barakuda)'

f_out.close()



print cf_nc+' created!!!'

