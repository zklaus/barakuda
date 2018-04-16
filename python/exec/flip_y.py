#!/usr/bin/env python

import sys
import os
#import numpy as nmp
from netCDF4 import Dataset
from string import replace

if len(sys.argv) != 2:
    print 'Usage: '+sys.argv[0]+' <runoff_depth.nc>'
    sys.exit(0)
cf_orig  = sys.argv[1]
#cf_in   = sys.argv[2]
#cmin_dept    = sys.argv[3] ; rmin_depth = float(cmin_dept)

cf_in = replace(cf_orig, '.nc', '_vfliped.nc')
os.system('rm -f '+cf_in)
os.system('cp '+cf_orig+' '+cf_in)

print '\n'


f_in = Dataset(cf_in, 'r+')

list_var = f_in.variables.keys()

nbvar = len(list_var)

print '\n *** '+str(nbvar)+'variables =>', list_var


for jv in range(nbvar):
    print ''
    
    cv = list_var[jv]
    l_dim = f_in.variables[cv].dimensions
    nbdim = len(l_dim)
    l_y = ( 'y' in l_dim )
    idx_y = -1
    if l_y: idx_y = l_dim.index('y')
    print '     '+cv+' => ', l_dim, nbdim, l_y, idx_y

    if l_y:

        if nbdim==4:
            xf = f_in.variables[cv][:,:,:,:]
            if idx_y==2:
                print '    flipping variable '+cv+' along axis # 2'
                f_in.variables[cv][:,:,:,:] = f_in.variables[cv][:,:,::-1,:]

        elif nbdim==3:
            xf = f_in.variables[cv][:,:,:]
            if idx_y==1:
                print '    flipping variable '+cv+' along axis # 1'
                f_in.variables[cv][:,:,:] = f_in.variables[cv][:,::-1,:]

        elif nbdim==2:
            xf = f_in.variables[cv][:,:]
            print '    flipping variable '+cv+' along axis # 0'
            if idx_y==0: f_in.variables[cv][:,:] = f_in.variables[cv][::-1,:]

        else:
            print 'ERROR: only dimensions 2 up to 4!' ; sys.exit(0)

        
f_in.close()

print '\n *** '+cf_in+' written!\n'
