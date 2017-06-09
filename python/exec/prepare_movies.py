#!/usr/bin/env python

#       B a r a K u d a
#
#  Prepare 2D maps (monthly) that will later become a GIF animation!
#  NEMO output and observations needed
#
#    L. Brodeau, november 2016

import sys
import os
import numpy as nmp
from string import replace
from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_plot as bp

venv_needed = {'ORCA','EXP','DIAG_D','MM_FILE','NN_SST','NN_T','NN_SSS','NN_S','NN_MLD',
               'FILE_ICE_SUFFIX','NN_ICEF',
               'NM_TS_OBS','F_T_OBS_3D_12','F_S_OBS_3D_12','F_SST_OBS_12','NN_SST_OBS','NN_T_OBS','NN_S_OBS'}
# 'NM_IC_OBS'

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFEXP = vdic['ORCA']+'-'+vdic['EXP']

tmin=-4.  ;  tmax=-tmin ;  dtemp = 0.25
smin=-1.4 ;  smax=-smin ;  dsali = 0.1
mmin=0.   ;  mmax=1500. ;  dmld  = 50.


fig_type='png'


cn_obs_ts = vdic['NM_TS_OBS']
#cn_obs_ic = vdic['NM_TS_OBS']

narg = len(sys.argv)
if narg < 4:
    print 'Usage: '+sys.argv[0]+' <NEMO file (1 year, monthyly)> <year> <var>'
    print '          with var one of "sst", "sss", "ice", "mld"'
    sys.exit(0)

cf_T = sys.argv[1]
cy    = sys.argv[2] ; jy=int(cy)
cvar  = sys.argv[3]

cf_ice = replace(cf_T, 'grid_T', vdic['FILE_ICE_SUFFIX'])

print ' *** file to read '+vdic['NN_ICEF']+' from: '+cf_ice+'\n'

if not cvar in ['sst','sss','ice','mld']:
    print 'ERROR (prepare_movies.py): variable '+cvar+' not supported yet!'
    sys.exit(0)

path_fig = 'movies'

os.system("mkdir -p "+path_fig)




# 3D climatology :
# ------------


# Salinity
if cvar == 'sss':
    bt.chck4f(vdic['F_S_OBS_3D_12'])
    id_clim = Dataset(vdic['F_S_OBS_3D_12'])
    Vclim  = id_clim.variables[vdic['NN_S_OBS']][:,0,:,:]; print '(has ',Vclim.shape[0],' time snapshots)\n'
    id_clim.close()

# 2D SST obs :
if cvar == 'sst':
    cv_sst_obs = vdic['NN_SST_OBS']
    bt.chck4f(vdic['F_SST_OBS_12'])
    id_clim_sst = Dataset(vdic['F_SST_OBS_12'])
    nb_dim = len(id_clim_sst.variables[cv_sst_obs].dimensions)
    if nb_dim == 3:
        Vclim  = id_clim_sst.variables[cv_sst_obs][:,:,:]; print '(has ',Vclim.shape[0],' time snapshots)\n'
    elif nb_dim == 4:
        Vclim  = id_clim_sst.variables[cv_sst_obs][:,0,:,:]; print '(has ',Vclim.shape[0],' time snapshots)\n'
    else:
        print 'ERROR (prepare_movies.py): shape of '+cv_sst_obs+' in '+vdic['F_SST_OBS_12']+' is problematic!'
        sys.exit(0)
    id_clim_sst.close()

# Sea-ice concentration :
# => no clim used!

if cvar in ['sss','sst']:  ( nmn , nj0 , ni0 ) = Vclim.shape



# Getting land-sea mask and coordinates:
bt.chck4f(vdic['MM_FILE'])
id_mask = Dataset(vdic['MM_FILE'])
xlon  = id_mask.variables['nav_lon'][:,:]
xlat  = id_mask.variables['nav_lat'][:,:]
imask = id_mask.variables['tmask'][0,0,:,:]
id_mask.close()



# Getting NEMO variables:
# -----------------------

cf_in = cf_T
if cvar == 'ice' and cf_ice != cf_T:
    cf_in = cf_ice

bt.chck4f(cf_in)

id_in = Dataset(cf_in)
#list_var = id_in.variables.keys()
if cvar == 'sst':
    if vdic['NN_SST'] == 'thetao' or vdic['NN_SST'] == 'votemper' :   #lolo:bad !!! should check shape!!!
        Vnemo = id_in.variables[vdic['NN_SST']][:,0,:,:]
    else:
        Vnemo = id_in.variables[vdic['NN_SST']][:,:,:]
    cv = 'dsst'

if cvar == 'sss':
    if vdic['NN_SSS'] == 'so' or vdic['NN_SSS'] == 'vosaline' :   #lolo:bad !!! should check shape!!!
        Vnemo = id_in.variables[vdic['NN_SSS']][:,0,:,:]
    else:
        Vnemo = id_in.variables[vdic['NN_SSS']][:,:,:]
    cv = 'dsss'

if cvar == 'ice':
    if vdic['NN_ICEF'] == 'X':
        print 'ERROR (prepare_movies.py): you set "X" (missing) as the name for ice concentration in your conf file!'; sys.exit(0)
    Vnemo = id_in.variables[vdic['NN_ICEF']][:,:,:]

if cvar == 'mld':
    if vdic['NN_MLD'] == 'X':
        print 'ERROR (prepare_movies.py): you set "X" (missing) as the name for MLD in your conf file!'; sys.exit(0)
    Vnemo = id_in.variables[vdic['NN_MLD']][:,:,:]
id_in.close()




[ nt, nj, ni ] = Vnemo.shape

if nt != 12:
    print 'ERROR (prepare_movies.py): we expect 12 montly records in NEMO grid_T file!'
    sys.exit(0)
    
if (cvar in ['sss','sst']) and (nj != nj0 or ni != ni0):
    print 'ERROR (prepare_movies.py): NEMO file and clim do no agree in shape for '+cvar+'!'
    print '       clim => '+str(ni0)+', '+str(nj0)+', ('+vdic['F_T_OBS_3D_12']+')'
    print '       NEMO => '+str(ni)+', '+str(nj)
    sys.exit(0)

if cvar in ['sss','sst','mld']:
    # Creating 1D long. and lat.:
    ji_lat0 = nmp.argmax(xlat[nj-1,:])
    vlon = nmp.zeros(ni) ; vlon[:] = xlon[20,:]
    vlat = nmp.zeros(nj) ; vlat[:] = xlat[:,ji_lat0]


if cvar == 'ice':
    # Extraoplating sea values over continents:
    bt.drown(Vnemo[:,:,:], imask, k_ew=2, nb_max_inc=10, nb_smooth=10)



lpix = False
if vdic['ORCA'][:5] == 'ORCA0': lpix = True


for jt in range(nt):

    cm = "%02d" % (jt+1)
    cdate  = cy+cm
    cdatet = cy+'/'+cm

    if cvar == 'sst':
        bp.plot("2d")(vlon, vlat, Vnemo[jt,:,:] - Vclim[jt,:,:],
                      imask[:,:],  tmin, tmax, dtemp,
                      corca=vdic['ORCA'], lkcont=False, cpal='RdBu_r',
                      cfignm=path_fig+'/'+cv+'_'+cdate,
                      cbunit='K', cfig_type=fig_type, lat_min=-77., lat_max=75.,
                      ctitle='SST (NEMO - "'+cn_obs_ts+'"), '+CONFEXP+' ('+cdatet+')',
                      lforce_lim=True, i_cb_subsamp=2, lpix=lpix)

    if cvar == 'sss':
        bp.plot("2d")(vlon, vlat, Vnemo[jt,:,:] - Vclim[jt,:,:],
                      imask[:,:],  smin, smax, dsali,
                      corca=vdic['ORCA'], lkcont=False, cpal='PiYG_r',
                      cfignm=path_fig+'/'+cv+'_'+cdate,
                      cbunit='PSU', cfig_type=fig_type, lat_min=-77., lat_max=75.,
                      ctitle='SSS (NEMO - "'+cn_obs_ts+'"), '+CONFEXP+' ('+cdatet+')',
                      lforce_lim=True, i_cb_subsamp=2, lpix=lpix)

    if cvar == 'mld':
        bp.plot("2d")(vlon, vlat, Vnemo[jt,:,:], imask[:,:],  mmin, mmax, dmld,
                      corca=vdic['ORCA'], lkcont=False, cpal='ncview_nrl',
                      cfignm=path_fig+'/'+cvar+'_'+cdate,
                      cbunit='m', cfig_type=fig_type, lat_min=-77., lat_max=75.,
                      ctitle='Mixed-Layer depth, '+CONFEXP+' ('+cdatet+')',
                      lforce_lim=True, i_cb_subsamp=2, lpix=lpix)
    
    if cvar == 'ice':
        # Extraoplating sea values on continents:
        bt.drown(Vnemo[jt,:,:], imask, k_ew=2, nb_max_inc=10, nb_smooth=10)
        # ICE north:
        cv = "icen"
        bp.plot("nproj")('npol2', 0., 1., 0.1, xlon, xlat, Vnemo[jt,:,:],
                         cfignm=path_fig+'/'+cv+'_'+cdate, cpal='ice', cbunit='(frac.)',
                         ctitle='Ice concentration, '+CONFEXP+' ('+cdatet+')',
                         lkcont=True, cfig_type=fig_type, lforce_lim=True)

        cv = "ices"
        bp.plot("nproj")('spstere', 0., 1., 0.1, xlon, xlat, Vnemo[jt,:,:],
                         cfignm=path_fig+'/'+cv+'_'+cdate, cpal='ice', cbunit='(frac.)',
                         ctitle='Ice concentration, '+CONFEXP+' ('+cdatet+')',
                         lkcont=True, cfig_type=fig_type, lforce_lim=True)





print '\n *** EXITING prepare_movies.py for year '+cy+', var ='+cvar+' !\n'
