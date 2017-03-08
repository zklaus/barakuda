#!/usr/bin/env python

#       B a r a K u d a
#
#       L. Brodeau, 2017]

import sys
import numpy as nmp
import string
import os
from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_orca as bo
import barakuda_ncio as bn

#cf_tf = '/home/laurent/DEV/barakuda/data/TS_sections.dat'

narg = len(sys.argv)
if narg != 3:
    print 'Usage: '+sys.argv[0]+' <TS_transect_file> <mesh_mask>'; sys.exit(0)
cf_tf = sys.argv[1]
cf_mm = sys.argv[2]

# Opening mesh_mask:
f_mm = Dataset(cf_mm)
Xlon = f_mm.variables['glamt'][0,:,:]
Xlat = f_mm.variables['gphit'][0,:,:]
mask    = f_mm.variables['tmask'][0,0,:,:]
f_mm.close()

(nj,ni) = nmp.shape(Xlon)

nmask = nmp.zeros((nj,ni))



# Getting sections:
vboxes, vlon1, vlat1, vlon2, vlat2 = bt.read_coor(cf_tf, ctype='float', lTS_bounds=False)


js = -1

for csname in vboxes:

    js = js + 1

    print'\n *** '+sys.argv[0]+': treating section '+csname


    nmask[:,:] = mask[:,:]

    ( i1, i2, j1, j2 ) = bo.transect_zon_or_med(vlon1[js], vlon2[js], vlat1[js], vlat2[js], Xlon, Xlat)


    print csname+' :'
    print   '(lon1, lon2, lat1, lat2) =', vlon1[js], vlon2[js], vlat1[js], vlat2[js]
    print   ' => i1, i2, j1, j2 =', i1, i2, j1, j2

    # Just for info:
    lat1_o = Xlat[j1,i1]
    lat2_o = Xlat[j2,i2]
    
    lon1_o = Xlon[j1,i1]
    lon2_o = Xlon[j2,i2]

    # Mid point:
    im = (i1+i2)/2
    jm = (j1+j2)/2
    lat_m = Xlat[jm,im]
    lon_m = Xlon[jm,im]

    if lon1_o > 180.: lon1_o = lon1_o - 360.
    if lon2_o > 180.: lon2_o = lon2_o - 360.
    if lon_m  > 180.: lon_m  = lon_m  - 360.

    print   ' => ORCA p1 , p2   =>' , lon1_o, lon2_o, lat1_o, lat2_o
    print   ' => ORCA mid point =>' , lon_m, lat_m
    print ''
    
    if i1 > i2: print 'ERROR: cross_sections.py => i1 > i2 !'; sys.exit(0)
    if j1 > j2: print 'ERROR: cross_sections.py => j1 > j2 !'; sys.exit(0)

    ip = 0
    if i1 == i2:
        print 'Meridional section!'
        ip = 1

    jp=0
    if j1 == j2:
        print 'Zonal section!'
        jp=1

    nmask[j1:j2+jp,i1:i2+ip] = -1
    
    #bn.write_2d_mask('mask_'+csname+'.nc', nmask, xlon=Xlon, xlat=Xlat)
    bn.write_2d_mask('mask_'+csname+'.nc', nmask)
        

