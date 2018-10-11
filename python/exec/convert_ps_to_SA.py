#!/usr/bin/env python


# L. Brodeau, 2015


# Practical salinity to absolute salinity (TEOS 10)

import sys
import os
import numpy as nmp
from netCDF4 import Dataset
from string import replace

import gsw

SSO = 35.16504

cdepth = 'deptht'

if len(sys.argv) != 3:
    print 'Usage: '+sys.argv[0]+' <Salinity_file_to_convert> <salinity_name>'
    sys.exit(0)


cf_sal  = sys.argv[1]
cv_sal  = sys.argv[2]

cf_out = replace(cf_sal, cv_sal, cv_sal+'_TEOS10')


os.system('rm -f '+cf_out)
os.system('cp '+cf_sal+' '+cf_out)


l_accurate = False




print '\n'

# Opening the Netcdf file:
f_out = Dataset(cf_out, 'r+')     # r+ => can read and write in the file... )
print 'File ', cf_out, 'is open...\n'

# Inquire variables in the file to see if a depth is there...
list_var = f_out.variables.keys()
print ' list_var =', list_var


if cdepth in list_var:
    l_accurate = True
    print ' *** Going for accurate method cause '+cdepth+' is present in '+cf_out+' !\n'


vcdim = f_out.variables[cv_sal].dimensions
cv_t = vcdim[0]; print ' *** record dimension is called "'+cv_t+'"'
Nt = f_out.dimensions[cv_t].size ; print ' *** There are '+str(Nt)+' time records...\n'

# Inquire the shape of arrays:
nb_dim = len(vcdim)
print ' *** '+cf_out+' has '+str(nb_dim)+' dimmensions!'

#if not nb_dim in [ 2, 3, 4 ]: print ' ERROR: unsupported number of dimmensions! =>', nb_dim ; sys.exit(0)
if not nb_dim in [ 4 ]: print ' ERROR: unsupported number of dimmensions! =>', nb_dim ; sys.exit(0)


if nb_dim != 4 and l_accurate:
    print '\n WARNING!!! => using less accurate method because nb_dim = '+str(nb_dim)
    print '          and we do not know how to handle this case yet....'
    l_accurate = False



for jt in range(Nt):
         
    print '\n --- treating record # '+str(jt)
    
    if nb_dim==4:
        xsal = f_out.variables[cv_sal][jt,:,:,:]
        if jt == 0:
            (nk,nj,ni) = nmp.shape(xsal)
        
    if nb_dim==3:
        xsal = f_out.variables[cv_sal][jt,:,:]
        if jt == 0: (nj,ni) = nmp.shape(xsal)
        
    if nb_dim==2:
        xsal = f_out.variables[cv_sal][jt,:]
        if jt == 0: (ni) = nmp.shape(xsal)

    if jt == 0: shp  = nmp.shape(xsal)


    xtmp = nmp.zeros(shp)

    xtmp = xsal
    
    if l_accurate:

        if jt == 0:
            print '\n Using accurate method with depth !'        
            vz    = f_out.variables['deptht'][:]
            xdepth = nmp.zeros(shp)
            # building 3d arrays of depth:
            for jk in range(nk):
                if nb_dim==4: xdepth[jk,:,:] = vz[jk]
                if nb_dim==3: xdepth[jk,:]   = vz[jk]
                if nb_dim==2: xdepth[jk]     = vz[jk]

        # pressure should be in dbar and it's the same as the depth in metre actually:
        if nb_dim==4: f_out.variables[cv_sal][jt,:,:,:] = gsw.SA_from_SP(xtmp[:,:,:], xdepth, -140., 0.)
        if nb_dim==3: f_out.variables[cv_sal][jt,:,:]   = gsw.SA_from_SP(xtmp[:,:],   xdepth, -140., 0.)
        if nb_dim==2: f_out.variables[cv_sal][jt,:]     = gsw.SA_from_SP(xtmp[:],     xdepth, -140., 0.)
    
    else:
    
        # Fabien says it's enough:
        if nb_dim==4: f_out.variables[cv_sal][jt,:,:,:] = xtmp[:,:,:]*SSO/35.
        if nb_dim==3: f_out.variables[cv_sal][jt,:,:]   = xtmp[:,:]*SSO/35.
        if nb_dim==2: f_out.variables[cv_sal][jt,:]     = xtmp[:]*SSO/35.


    
        
f_out.variables[cv_sal].long_name = 'Absolute Salinity (TEOS10) build from practical salinity'
    
f_out.close()

print cf_out+' sucessfully created!'

