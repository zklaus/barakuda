#!/usr/bin/env python
#
# L. Brodeau, 2017

#
# Compute curl, aka relative vorticity from 2D vector on a lat-lon spherical
# grid a la ECMWF
#

import sys
import numpy as nmp
from netCDF4 import Dataset
import string

import barakuda_tool as bt


cv_Vx='EWSS' ;# cv_Vx_i='Taux'
cv_Vy='NSSS' ;# cv_Vy_i='Tauy'
cv_tm='CURL'

radius_Earth = 6371. ; #  km !!! Not m !!! => curl will be in XXX !!!


to_rad = nmp.pi/180.


#jt1=10
#jt2=12

if len(sys.argv) != 2:
    print 'Usage: '+sys.argv[0]+' <FILE_TAU_X.nc>'
    sys.exit(0)

cf_Vx   = sys.argv[1]
cf_Vy = string.replace(cf_Vx, cv_Vx, cv_Vy)
cf_crl = string.replace(cf_Vx, cv_Vx, cv_tm)


print " cf_Vx = ", cf_Vx
print " cf_Vy = ", cf_Vy
print " cf_crl = ", cf_crl
print "\n"


#  TAUY
#  ~~~~
bt.chck4f(cf_Vy) ; f_Vy_in = Dataset(cf_Vy)

# Extracting the longitude and 1D array:
vlon     = f_Vy_in.variables['lon'][:]
clnm_lon = f_Vy_in.variables['lon'].long_name ; cunt_lon = f_Vy_in.variables['lon'].units
csnm_lon = f_Vy_in.variables['lon'].standard_name
print 'LONGITUDE: ', clnm_lon, cunt_lon, csnm_lon

# Extracting the longitude 1D array:
vlat     = f_Vy_in.variables['lat'][:]
clnm_lat = f_Vy_in.variables['lat'].long_name ; cunt_lat = f_Vy_in.variables['lat'].units
csnm_lat = f_Vy_in.variables['lat'].standard_name
print 'LATGITUDE: ', clnm_lat, cunt_lat, csnm_lat

# Extracting time 1D array:
vtime     = f_Vy_in.variables['time'][:] ; cunt_time = f_Vy_in.variables['time'].units
print 'TIME: ', cunt_time, '\n'

# Extracting a variable, ex: "t" the 3D+T field of temperature:
xty     = f_Vy_in.variables[cv_Vy][:,:,:]
cunt_Vy = f_Vy_in.variables[cv_Vy].units
code_Vy = f_Vy_in.variables[cv_Vy].code
ctab_Vy = f_Vy_in.variables[cv_Vy].table
print cv_Vy+': ', cunt_Vy, code_Vy, ctab_Vy, '\n'
f_Vy_in.close()



# TAUX
# ~~~
bt.chck4f(cf_Vx) ; f_Vx_in = Dataset(cf_Vx)
xtx     = f_Vx_in.variables[cv_Vx][:,:,:]
cunt_Vx = f_Vx_in.variables[cv_Vx].units
code_Vx = f_Vx_in.variables[cv_Vx].code
ctab_Vx = f_Vx_in.variables[cv_Vx].table
print cv_Vx+': ', cunt_Vx, code_Vx, ctab_Vx, '\n'
f_Vx_in.close()



# Checking dimensions
# ~~~~~~~~~~~~~~~~~~~
dim_Vy = xty.shape ; dim_Vx = xtx.shape
if dim_Vy != dim_Vx:
    print 'Shape problem!!!'; print dim_Vy , dim_Vx

print '\n'
( Nt, nj, ni ) = dim_Vy
print 'ni, nj, Nt = ', ni, nj, Nt



# Building tm
# ~~~~~~~~~~~

xcurl = nmp.zeros(( Nt, nj, ni ))
#xcurl = nmp.sqrt( xtx*xtx + xty*xty )



dlamx2 = 2.*(vlon[1] - vlon[0])*to_rad
dphix2 = 2.*(vlat[0] - vlat[1])*to_rad

xcosphi = nmp.zeros(( nj, ni ))
one_on_xcosphi = nmp.zeros(( nj, ni ))
for ji in range(ni):
    xcosphi[:,ji] = nmp.cos(vlat[:]*to_rad)

one_on_xcosphi[:,:] = 1./xcosphi[:,:]

print ' dlamx2, dphix2 = ', dlamx2, dphix2


xtmp = nmp.zeros(( nj, ni ))


# Creating output file
# ~~~~~~~~~~~~~~~~~~~~
f_out = Dataset(cf_crl, 'w', format='NETCDF3_CLASSIC')

# Dimensions:
f_out.createDimension('lon', ni)
f_out.createDimension('lat', nj)
f_out.createDimension('time', None)

# Variables
id_lon = f_out.createVariable('lon','f4',('lon',))
id_lat = f_out.createVariable('lat','f4',('lat',))
id_tim = f_out.createVariable('time','f4',('time',))
id_tm  = f_out.createVariable(cv_tm,'f4',('time','lat','lon',))

# Attributes
id_tim.units = cunt_time

id_lat.long_name     = clnm_lat
id_lat.units         = cunt_lat
id_lat.standard_name = csnm_lat

id_lon.long_name     = clnm_lon
id_lon.units         = cunt_lon
id_lon.standard_name = csnm_lon

id_tim.units         = cunt_time

id_tm.long_name = 'Module of surface wind stress TAUY and TAUX'
id_tm.units = 'N/m^2'
id_tm.code  = '???'
id_tm.table = '128'

f_out.About = 'Created by L. Brodeau using TAUY and TAUX corresponding fields'

# Filling variables:
id_lat[:] = vlat[:]
id_lon[:] = vlon[:]




for jt in range(Nt):


    print ' *** jt = ', jt

    xtmp[:,:] = xtx[jt,:,:] * xcosphi[:,:]

    xcurl[jt,1:nj-1,1:ni-1] =   ( xty[jt,1:nj-1,2:ni]   - xty[jt,1:nj-1,0:ni-2] ) / dlamx2 \
    - one_on_xcosphi[1:nj-1,1:ni-1]*( xtmp[2:nj,1:ni-1] -   xtmp[0:nj-2,1:ni-1] ) / dphix2

    xcurl[jt,:,:] = -1./radius_Earth * xcurl[jt,:,:]
                                                
    id_tim[jt]     = vtime[jt]
    id_tm[jt,:,:]  = xcurl[jt,:,:] 

f_out.close()

print 'Bye!'

