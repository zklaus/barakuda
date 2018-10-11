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


# Inquire the shape of arrays:
nb_dim = len(f_out.variables[cv_sal].dimensions)
print ' *** '+cf_out+' has '+str(nb_dim)+' dimmensions!'
if   nb_dim==4:
    xsal = f_out.variables[cv_sal][:,:,:,:]
elif nb_dim==3:
    xsal = f_out.variables[cv_sal][:,:,:]
elif nb_dim==2:
    xsal = f_out.variables[cv_sal][:,:]
else:
    print ' ERROR: unsupported number of dimmensions!' ; sys.exit(0)




if nb_dim != 4 and l_accurate:
    print '\n WARNING!!! => using less accurate method because nb_dim = '+str(nb_dim)
    print '          and we do not know how to handle this case yet....'
    l_accurate = False
    

if l_accurate:

    print '\n Using accurate method with depth !'
    
    vz    = f_out.variables['deptht'][:]

    [nt,nk,nj,ni] = nmp.shape(xsal)

    print ' [nt,nk,nj,ni] =>', nt,nk,nj,ni
    
    xdepth = nmp.zeros((nk,nj,ni))

    # building 3d arrays of depth:
    for jk in range(nk): xdepth[jk,:,:] = vz[jk]

    # pressure should be in dbar and it's the same as the depth in metre actually:
    for jt in range(nt):
        print ' jt =', jt
        f_out.variables[cv_sal][jt,:,:,:] = gsw.SA_from_SP(xsal[jt,:,:,:], xdepth, -140., 0.)

else:

    # Fabien says it's enough:
    if nb_dim==2:
        f_out.variables[cv_sal][:,:]   = xsal[:,:]*SSO/35.
    if nb_dim==3:
        f_out.variables[cv_sal][:,:,:]   = xsal[:,:,:]*SSO/35.
    elif nb_dim==4:
        f_out.variables[cv_sal][:,:,:,:] = xsal[:,:,:,:]*SSO/35.
    else:
        print ' WTF???'; sys.exit(1)


f_out.variables[cv_sal].long_name = 'Absolute Salinity (TEOS10) build from practical salinity'

f_out.close()




print cf_out+' sucessfully created!'

