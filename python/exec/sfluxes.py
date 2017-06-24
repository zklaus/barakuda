#!/usr/bin/env python

#       B a r a K u d a
#
#     Generate zonal (and maybe one day global) plots related to
#     air-sea fluxes
#     Use climatology fields built with 'build_clim.sh'
#
#       L. Brodeau, 2017

import sys
from os import path
import numpy as nmp
from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_plot as bp

  
venv_needed = {'ORCA','EXP','DIAG_D','MM_FILE','FIG_FORM','FILE_FLX_SUFFIX',
               'NN_QNET','NN_QSOL','NM_QSOL_OBS','F_QSOL_OBS_12','NN_QSOL_OBS'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFEXP = vdic['ORCA']+'-'+vdic['EXP']

cd_clim = vdic['DIAG_D']+'/clim'

path_fig='./'
fig_type=vdic['FIG_FORM']

narg = len(sys.argv)
if narg < 3: print 'Usage: '+sys.argv[0]+' <year1> <year2>'; sys.exit(0)
cy1 = sys.argv[1] ; cy2=sys.argv[2]; jy1=int(cy1); jy2=int(cy2)
jy1_clim = jy1 ; jy2_clim = jy2
print ' => mean on the clim : ', jy1_clim, jy2_clim, '\n'

ctag = CONFEXP+'_'+cy1+'-'+cy2

ct_nemo = vdic['FILE_FLX_SUFFIX'] ; # (SBC or grid_T, ect.)


# Getting coordinates and mask:
bt.chck4f(vdic['MM_FILE'])
id_mm = Dataset(vdic['MM_FILE'])
xlon   = id_mm.variables['glamt'][0,:,:]
xlat = id_mm.variables['gphit'][0,:,:]
Xmask = id_mm.variables['tmask'][0,0,:,:]
id_mm.close()



#  Getting NEMO mean monthly climatology surface flux:
cextra_qsol = ''
cv_qsol = vdic['NN_QSOL']
if cv_qsol != 'X':
    cf_nemo_qsol = cd_clim+'/mclim_'+ctag+'_'+ct_nemo+'.nc4'
else:
    print '\n *** WARNING: sfluxes.py : if no "NN_QSOL" ! No sfluxes.py! Aborting !!!\n'
    sys.exit(0)
    
bt.chck4f(cf_nemo_qsol)
id_nemo = Dataset(cf_nemo_qsol)
Xqsol_nemo   = id_nemo.variables[cv_qsol][:,:,:]
id_nemo.close()

[ Nt, nj, ni ] = Xqsol_nemo.shape ; print ' Shape of QSOL :', Nt, nj, ni, '\n'
if not Nt in [1,12]:
    print '\n *** ERROR: sfluxes.py => only accepting monthly or annual climatologies (NEMO)!', Nt
    sys.exit(0)
if ( nj, ni ) != nmp.shape(Xmask):
    print '\n *** ERROR: sfluxes.py => mask and field disagree in shape!'; sys.exit(0)





# Getting OBS.:
cf_obs_qsol = vdic['F_QSOL_OBS_12'] ; cv_obs_qsol = vdic['NN_QSOL_OBS']
bt.chck4f(cf_obs_qsol)
id_obs = Dataset(cf_obs_qsol)
Xqsol_obs = id_obs.variables[cv_obs_qsol][:,:,:]
vlat_obs  = id_obs.variables['lat'][:]
id_obs.close()

[ Nto, njo, nio ] = Xqsol_obs.shape ; print ' Shape of observation QSOL :', Nto, njo, nio, '\n'
if not Nto == Nt:
    print '\n *** ERROR: sfluxes.py => Obs. and NEMO clim should have same number of time records!', Nto, Nt
    sys.exit(0)



# Zonal stuffs....
vz_obs = nmp.zeros((Nt,njo)) ; # a zonal profile...
vz_obs[:,:] = bt.mk_zonal(Xqsol_obs, r_mask_from_val=-9999.)


# Creating 1D long. and lat. for NEMO
vlon = nmp.zeros(ni) ; vlon[:] = xlon[nj/8,:]
ji_lat0 = nmp.argmax(xlat[nj-1,:])
vlat = nmp.zeros(nj) ; vlat[:] = xlat[:,ji_lat0]

vz_nemo = nmp.zeros((Nt,nj)) ; # a zonal profile...
vz_nemo[:,:] = bt.mk_zonal(Xqsol_nemo, XMSK=Xmask)

ctt = 'Zonnaly-averaged annual Net Solar heat flux, '+CONFEXP+' ('+cy1+'-'+cy2+') vs OBS.'
bp.plot("zonal")(vlat, nmp.mean(vz_nemo,axis=0), VY1=vlat_obs, VZ1=nmp.mean(vz_obs,axis=0),
                 cfignm=path_fig+'zonal_Qsol_vs_obs_annual_'+CONFEXP, zmin=50., zmax=250., dz=10.,
                 xmin=-70., xmax=65., dx=10., czunit=r'$(W/m^2)$', cfig_type=fig_type,
                 ctitle=ctt, lab='NEMO ('+cv_qsol+')', lab1=vdic['NM_QSOL_OBS'], loc_legend='center')


if Nt == 12:

    ctt = 'Zonnaly-averaged JFM Net Solar heat flux, '+CONFEXP+' ('+cy1+'-'+cy2+') vs OBS.'
    bp.plot("zonal")(vlat, nmp.mean(vz_nemo[:3,:],axis=0), VY1=vlat_obs, VZ1=nmp.mean(vz_obs[:3,:],axis=0),
                     cfignm=path_fig+'zonal_Qsol_vs_obs_JFM_'+CONFEXP, zmin=0., zmax=260., dz=20.,
                     xmin=-70., xmax=65., dx=10., czunit=r'$(W/m^2)$', cfig_type=fig_type,
                     ctitle=ctt, lab='NEMO ('+cv_qsol+')', lab1=vdic['NM_QSOL_OBS'], loc_legend='center')

    ctt = 'Zonnaly-averaged JAS Net Solar heat flux, '+CONFEXP+' ('+cy1+'-'+cy2+') vs OBS.'
    bp.plot("zonal")(vlat, nmp.mean(vz_nemo[6:9,:],axis=0), VY1=vlat_obs, VZ1=nmp.mean(vz_obs[6:9,:],axis=0),
                     cfignm=path_fig+'zonal_Qsol_vs_obs_JAS_'+CONFEXP, zmin=0., zmax=260., dz=20.,
                     xmin=-70., xmax=65., dx=10., czunit=r'$(W/m^2)$', cfig_type=fig_type,
                     ctitle=ctt, lab='NEMO ('+cv_qsol+')', lab1=vdic['NM_QSOL_OBS'], loc_legend='center')
