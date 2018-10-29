#!/usr/bin/env python


# L. Brodeau, 2015


# Potential temperature to conservative temperature (TEOS 10)

import sys
import os
import numpy as nmp
from netCDF4 import Dataset
from string import replace

import gsw

if len(sys.argv) != 5:
    print 'Usage: '+sys.argv[0]+' <Temperature_file_to_convert> <temperature_name> <Absolute_salinity_file> <salinity_name>'
    sys.exit(0)


cf_temp  = sys.argv[1]
cv_temp  = sys.argv[2]
cf_sal   = sys.argv[3]
cv_sal   = sys.argv[4]

cf_out = replace(cf_temp, cv_temp, cv_temp+'_TEOS10')

os.system('rm -f '+cf_out)
os.system('cp '+cf_temp+' '+cf_out)







print '\n'


f_sal = Dataset(cf_sal)     # r+ => can read and write in the file... )


vcdim = f_sal.variables[cv_sal].dimensions
cv_t = vcdim[0]; print ' *** record dimension is called "'+cv_t+'"'
Nt = f_sal.dimensions[cv_t].size ; print ' *** There are '+str(Nt)+' time records...\n'

# Inquire the shape of arrays:
nb_dim = len(vcdim)
print ' *** '+cf_sal+' has '+str(nb_dim)+' dimmensions!'

if not nb_dim in [ 2, 3, 4 ]: print ' ERROR: unsupported number of dimmensions! =>', nb_dim ; sys.exit(0)


# Opening the Netcdf output file:
f_out = Dataset(cf_out, 'r+')     # r+ => can read and write in the file... )
print 'File ', cf_out, 'is open...\n'



for jt in range(Nt):

    print '\n --- treating record # '+str(jt)

    if nb_dim==4: xsal = f_sal.variables[cv_sal][jt,:,:,:]
    if nb_dim==3: xsal = f_sal.variables[cv_sal][jt,:,:]
    if nb_dim==2: xsal = f_sal.variables[cv_sal][jt,:]

    # Extracting tmask at surface level:
    if nb_dim==4:
        xtemp = f_out.variables[cv_temp][jt,:,:,:]
        f_out.variables[cv_temp][jt,:,:,:] = gsw.CT_from_pt(xsal, xtemp)

    if nb_dim==3:
        xtemp = f_out.variables[cv_temp][jt,:,:]
        f_out.variables[cv_temp][jt,:,:] = gsw.CT_from_pt(xsal, xtemp)

    if nb_dim==2:
        xtemp = f_out.variables[cv_temp][jt,:]
        f_out.variables[cv_temp][jt,:] = gsw.CT_from_pt(xsal, xtemp)

    


f_out.variables[cv_temp].long_name = 'Conservative Temperature (TEOS10) built from potential temperature'

f_sal.close()
f_out.close()



print cf_out+' sucessfully created!'

