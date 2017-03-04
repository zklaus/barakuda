#!/usr/bin/env python
#
#       B a r a K u d a
#
#     Generate misc. spatial 3D averaging out of NEMO output files...
#
#       L. Brodeau, november 2013
#

import sys
import os
import numpy as nmp

from netCDF4 import Dataset
from string  import replace

import barakuda_tool as bt
import barakuda_orca as bo
import barakuda_ncio as bnc


venv_needed = {'ORCA','EXP','DIAG_D','MM_FILE','BM_FILE',
               'NEMO_SAVED_FILES','NN_T','NN_S','ANNUAL_3D','TSTAMP'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFEXP = vdic['ORCA']+'-'+vdic['EXP']

if len(sys.argv) != 4:
    print 'Usage : sys.argv[1] <ORCA1_EXP_grid_T.nc> <year> <T or S>'; sys.exit(0)

cnexec  = sys.argv[0]
cf_T_in = sys.argv[1]
cyear   = sys.argv[2] ; jyear = int(cyear); cyear = '%4.4i'%jyear
cv_todo = sys.argv[3]

print 'Current year is '+cyear+' !\n'


if not cv_todo in ['T','S']:
    print 'Usage : sys.argv[1] <ORCA1_EXP_grid_T.nc> <year> <T or S>'; sys.exit(0)
if cv_todo == 'T': cvar = vdic['NN_T']
if cv_todo == 'S': cvar = vdic['NN_S']


# Checking if the land-sea mask file is here:
for cf in [vdic['MM_FILE'], vdic['BM_FILE']]:
    if not os.path.exists(cf):
        print 'Mask file '+cf+' not found'; sys.exit(0)

# Reading the grid metrics:
id_mm = Dataset(vdic['MM_FILE'])
list_variables = id_mm.variables.keys()
rmask  = id_mm.variables['tmask'][0,:,:,:]
xe1t   = id_mm.variables['e1t'][0,:,:]
xe2t   = id_mm.variables['e2t'][0,:,:]
if 'e3t_0' in list_variables[:]:
    Xe3t = id_mm.variables['e3t_0'][0,:,:,:] # we need the 3D field becaus partial steps!!!
else:
    print 'ERROR: '+cnexec+' => how do we retrieve 3D e3t???'; sys.exit(0)
id_mm.close()

[ nk, nj, ni ] = rmask.shape




Xa   = nmp.zeros((nj, ni))
Xv   = nmp.zeros((nk, nj, ni))
Xa[:,:] = xe1t[:,:]*xe2t[:,:]
del xe1t, xe2t

for jk in range(nk): Xv[jk,:,:] = Xa[:,:]*Xe3t[jk,:,:]
del Xe3t


print 'Opening different basin masks in file '+vdic['BM_FILE']
list_basin_names, list_basin_lgnms = bo.get_basin_info(vdic['BM_FILE'])
nb_basins = len(list_basin_names)
mask      = nmp.zeros((nb_basins,nk,nj,ni))
msk_tmp   = nmp.zeros((nj,ni))
mask[0,:,:,:] = rmask[:,:,:] ; # global

id_bm = Dataset(vdic['BM_FILE'])
for jb in range(1,nb_basins) :
    msk_tmp[:,:] = id_bm.variables['tmask'+list_basin_names[jb]][:,:]
    for jk in range(nk):
        mask[jb,jk,:,:] = msk_tmp[:,:]*rmask[jk,:,:]
id_bm.close()

del rmask, msk_tmp


#############################################
# 3D averaging for temperature and salinity #
#############################################

print '\n\n +++ '+cnexec+' => Starting 3D-averaging diags!'

if vdic['ANNUAL_3D'] == '1y':
    cf_T_in = replace(cf_T_in, vdic['TSTAMP'], vdic['ANNUAL_3D'])
    print '  ==> USING '+vdic['ANNUAL_3D']+' file !!! =>', cf_T_in

print '      ==> variable '+cvar

# DATA:
id_in = Dataset(cf_T_in)
vdepth = id_in.variables['deptht'][:]
Xd_m = id_in.variables[cvar][:,:,:,:]
id_in.close()

print '      ==> variable '+cvar+' read !'

j100m  = bt.find_index_from_value(100.  , vdepth) ; print 'j100m  = ', j100m,  '=> ', vdepth[j100m]
j1000m = bt.find_index_from_value(1000. , vdepth) ; print 'j1000m = ', j1000m, '=> ', vdepth[j1000m]

[ nt, nk0, nj0, ni0 ] = Xd_m.shape

if nt != 12 and nt != 1 : print 'ERROR: '+cnexec+' => only treating monthly or annual data so far...'; sys.exit(0)

if [ nk0, nj0, ni0 ] != [ nk, nj, ni ]: print 'ERROR: '+cnexec+' => Field and metrics do not agree in size!'; sys.exit(0)

vtime = nmp.zeros(nt)
for jt in range(nt): vtime[jt] = float(jyear) + (float(jt)+0.5)*1./float(nt)
print ' * Calendar: ', vtime[:]

# Annual mean array for current year:
Xd_y = nmp.zeros((1, nk, nj, ni))

if nt == 12:
    Xd_y[0,:,:,:] = nmp.mean(Xd_m, axis=0)
else:
    Xd_y[0,:,:,:] = Xd_m[0,:,:,:]

joce = 0

for cocean in list_basin_names[:]:

    colnm = list_basin_lgnms[joce]

    print '   ===> '+cvar+' in basin '+cocean+' ('+colnm+')'

    # Decrasing the domain size if possible:
    (vjj , vji)  = nmp.where(mask[joce,0,:,:]>0.5)
    j1 = max( nmp.min(vjj)-2 , 0    )
    i1 = max( nmp.min(vji)-2 , 0    )
    j2 = min( nmp.max(vjj)+2 , nj-1 ) + 1
    i2 = min( nmp.max(vji)+2 , ni-1 ) + 1

    if (i1,j1,i2,j2) != (0,0,ni,nj): print '       ===> zooming on i1,j1 -> i2,j2:', i1,j1,'->',i2,j2
    
    # I) Montly mean for diffrent depth ranges
    # ========================================
    
    Vts_tot      = bo.mean_3d(Xd_m[:,:,j1:j2,i1:i2],            mask[joce,:,j1:j2,i1:i2],            Xv[:,j1:j2,i1:i2]) ; # Top to bottom
    Vts_0_100    = bo.mean_3d(Xd_m[:,:j100m,j1:j2,i1:i2],       mask[joce,:j100m,j1:j2,i1:i2],       Xv[:j100m,j1:j2,i1:i2])
    Vts_100_1000 = bo.mean_3d(Xd_m[:,j100m:j1000m,j1:j2,i1:i2], mask[joce,j100m:j1000m,j1:j2,i1:i2], Xv[j100m:j1000m,j1:j2,i1:i2])
    Vts_1000_bot = bo.mean_3d(Xd_m[:,j1000m:,j1:j2,i1:i2],      mask[joce,j1000m:,j1:j2,i1:i2],      Xv[j1000m:,j1:j2,i1:i2])

    cf_out = vdic['DIAG_D']+'/3d_'+cvar+'_'+CONFEXP+'_'+cocean+'.nc'
    cv1 = cvar+'_0-bottom'
    cv2 = cvar+'_0-100'
    cv3 = cvar+'_100-1000'
    cv4 = cvar+'_1000-bottom'

    bnc.wrt_appnd_1d_series(vtime, Vts_tot, cf_out, cv1,
                            cu_t='year', cu_d='Unknown', cln_d ='3D-average of '+cvar+': surface to bottom, '+colnm,
                            vd2=Vts_0_100,    cvar2=cv2, cln_d2='3D-average of '+cvar+': surface to 100m, '+colnm,
                            vd3=Vts_100_1000, cvar3=cv3, cln_d3='3D-average of '+cvar+': 100m to 1000m, '+colnm,
                            vd4=Vts_1000_bot, cvar4=cv4, cln_d4='3D-average of '+cvar+': 1000m to bottom, '+colnm)



    # II) Annual mean vertical profile
    # ================================

    Vf = nmp.zeros(nk)

    for jk in range(nk):

        [ rf ] = bo.mean_2d(Xd_y[:,jk,j1:j2,i1:i2], mask[joce,jk,j1:j2,i1:i2], Xa[j1:j2,i1:i2])

        Vf[jk] = rf



    # NETCDF:
    cf_out = vdic['DIAG_D']+'/'+cvar+'_mean_Vprofile_'+CONFEXP+'_'+cocean+'.nc'
    l_nc_is_new = not os.path.exists(cf_out)
    #
    # Creating/Opening output Netcdf file:
    if l_nc_is_new:
        f_out = Dataset(cf_out, 'w', format='NETCDF3_CLASSIC')
    else:
        f_out = Dataset(cf_out, 'a', format='NETCDF3_CLASSIC')

    if l_nc_is_new:
        jrec2write = 0

        # Creating Dimensions:
        f_out.createDimension('time', None)
        f_out.createDimension('deptht', nk)

        # Creating variables:
        id_t = f_out.createVariable('time','f4',('time',)) ;      id_t.units = 'year'
        id_z = f_out.createVariable('deptht','f4',('deptht',)) ;  id_z.units = 'm'
        id_v01   = f_out.createVariable(cvar ,'f4',('time','deptht',))
        id_v01.long_name = 'Horizontally-averaged '+cvar+': '+colnm
        # Writing depth vector
        id_z[:] = vdepth[:]
        id_t[jrec2write] = float(jyear)+0.5
        id_v01[jrec2write,:] = Vf[:]
        f_out.Author = 'L. Brodeau ('+cnexec+' of Barakuda)'

    else:
        vt = f_out.variables['time']
        jrec2write = len(vt)
        v01 = f_out.variables[cvar]
        vt[jrec2write] = float(jyear)+0.5
        v01[jrec2write,:] = Vf[:]

    f_out.close()
    print cf_out+' written!'





    print ''

    joce = joce + 1


del Xd_m

print '\n'

print ' +++ '+cnexec+' => Done with 3D-averaging of variable '+cvar+'!\n'


print '\n *** EXITING '+cnexec+' for year '+cyear+' and variable '+cvar+'!\n'
