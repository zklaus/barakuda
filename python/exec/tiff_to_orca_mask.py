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

vbasins = [   'pac'   ,    'atl' ,     'ind'  ,   'soc'   ,   'arc'  ,  'wed'    ,  'lab'     ,  'med'  , 'gin'  ]
vbnames = [ 'Pacific' , 'Atlantic' , 'Indian' , 'Southern', 'Arctic' , 'Weddell' , 'Labrador' ,  'Med'  , 'GIN'  ]
vocesea = [ 'Ocean'   ,  'Ocean'   ,  'Ocean' ,  'Ocean'  ,  'Ocean' ,   'Sea'   ,   'Sea'    ,  'Sea'  , 'Seas' ]
vmandat = [  True     ,   True     ,   True   ,    False  ,    False ,   False   ,   False    ,   False ,  False ] ; # Mandatory ?

# Opening mesh_mask:
f_mm = Dataset(cf_mm)
nav_lon = f_mm.variables['nav_lon'][:,:]
nav_lat = f_mm.variables['nav_lat'][:,:]
mask    = f_mm.variables['tmask'][0,0,:,:]
f_mm.close()

(nj,ni) = nmp.shape(nav_lon)

cf_bm = string.replace(os.path.basename(cf_mm), 'mesh_', 'basin_')

nb_bas = len(vbasins)

vtreat = nmp.zeros(nb_bas, dtype=bool)

jb = 0 ; Nbt = 0
for cb in vbasins:

    cf_tiff = 'mask_'+cb+'.tiff'

    if not os.path.exists(cf_tiff):
        if vmandat[jb]:
            print 'PROBLEM! image '+cf_tiff+' is not here!!!' ; sys.exit(0)
    else:
        print ' *** good, '+cf_tiff+' is here...'
        vtreat[jb] = True
        Nbt = Nbt + 1

    jb = jb+1

#print vtreat


XBASINS = nmp.zeros((Nbt,nj,ni))

jb = 0 ; jbt = 0
for cb in vbasins:

    if vtreat[jb]:

        cf_tiff = 'mask_'+cb+'.tiff'

        # Opening Images:
        pic = Image.open(cf_tiff)

        xtmp = nmp.flipud( nmp.array(pic) )
        (njb,nib) = xtmp.shape

        if (njb,nib) != (nj,ni):
            print 'ERROR: something is wrong with the shapes of 2D arrays with basin '+cb+':'
            print nj,ni
            print njb,nib
            sys.exit(0)

        xmsk    = nmp.zeros((nj,ni))
        idx_sea = nmp.where(xtmp > 0)
        xmsk[idx_sea] = 1
        XBASINS[jbt,:,:] = xmsk[:,:]
        jbt = jbt + 1

    jb = jb + 1



now   = datetime.datetime.now()
cdate = now.strftime("%Y-%m-%d")



f_out = Dataset(cf_bm, 'w', format='NETCDF3_CLASSIC')

# Dimensions:
f_out.createDimension('x', ni)
f_out.createDimension('y', nj)


id_lon  = f_out.createVariable('nav_lon' ,'f4',('y','x',))
id_lat  = f_out.createVariable('nav_lat' ,'f4',('y','x',))

id_lon[:,:] = nav_lon[:,:]
id_lat[:,:] = nav_lat[:,:]


jb = 0 ; jbt = 0
for cb in vbasins:

    if vtreat[jb]:
        id_bas  = f_out.createVariable('tmask'+cb,'f4',('y','x',))
        id_bas.long_name = vbnames[jb]+' '+vocesea[jb]+' basin'
        id_bas[:,:] = XBASINS[jbt,:,:]*mask[:,:]
        jbt = jbt + 1
        
    jb = jb + 1

f_out.About  = 'ORCA025, masks for main ocean basins, created with orca_mesh_mask_to_bitmap.py, Gimp, and tiff_to_orca_mask.py, '+cdate+'.'
f_out.Author = 'L. Brodeau (https://github.com/brodeau/barakuda)'

f_out.close()



print cf_bm+' created!!!'

