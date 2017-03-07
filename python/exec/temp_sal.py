#!/usr/bin/env python

#       B a r a K u d a
#
#  Generate misc. spatial plots of potential temperature and salinity out of
#  NEMO output and climatology (from initial condition and surface restoring
#  NEMO files)
#
#    L. Brodeau, november 2013

import sys
import os
import numpy as nmp

from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_orca as bo
import barakuda_plot as bp

lfig0 = True
lfig1 = True
#lfig0 = False
#lfig1 = False
lfig2 = True

venv_needed = {'ORCA','EXP','DIAG_D','COMP2D','i_do_sect','MM_FILE','NN_SST','NN_T','NN_SSS','NN_S',
               'F_T_CLIM_3D_12','F_S_CLIM_3D_12','SST_CLIM_12','NN_SST_CLIM','NN_T_CLIM','NN_S_CLIM'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFEXP = vdic['ORCA']+'-'+vdic['EXP']

CC = vdic['COMP2D']

i_do_sect = int(vdic['i_do_sect'])


# Bounds and increment for comparison maps:
if CC == 'CLIM':
    tmin=-5.  ;  tmax=5.  ; dtemp = 0.5
    smin=-2.5 ;  smax=2.5 ; dsali = 0.25
else:
    tmin=-1.  ;  tmax=1. ;  dtemp = 0.05
    smin=-0.5 ;  smax=.5 ;  dsali = 0.025


path_fig='./'


fig_type='png'



narg = len(sys.argv)
if narg < 3: print 'Usage: '+sys.argv[0]+' <year1> <year2>'; sys.exit(0)
cy1 = sys.argv[1] ; cy2=sys.argv[2]; jy1=int(cy1); jy2=int(cy2)

#

if not ( jy1 >= 1984 and jy2 <= 2006 ):
    jy1_clim = 1984 ; jy2_clim = 2006
else:
    jy1_clim = jy1 ;  jy2_clim = jy2

print ' First and last year to treat:', jy1, jy2
print ' => mean on the clim : ', jy1_clim, jy2_clim, '\n'




# 3D climatology :
# ------------

# Temperature
bt.chck4f(vdic['F_T_CLIM_3D_12'])
id_clim = Dataset(vdic['F_T_CLIM_3D_12'])
Tclim  = id_clim.variables[vdic['NN_T_CLIM']][:,:,:,:]; print '(has ',Tclim.shape[0],' time snapshots)\n'
id_clim.close()
[ nmn , nk0 , nj0 , ni0 ] = Tclim.shape

# Salinity
bt.chck4f(vdic['F_S_CLIM_3D_12'])
id_clim = Dataset(vdic['F_S_CLIM_3D_12'])
Sclim  = id_clim.variables[vdic['NN_S_CLIM']][:,:,:,:]; print '(has ',Sclim.shape[0],' time snapshots)\n'
id_clim.close()




# 2D SST obs :
print 'We use the following SST climatology:'; print vdic['SST_CLIM_12']
bt.chck4f(vdic['SST_CLIM_12'])
id_clim_sst = Dataset(vdic['SST_CLIM_12'])
SSTclim  = id_clim_sst.variables[vdic['NN_SST_CLIM']][:,:,:]; print '(has ',SSTclim.shape[0],' time snapshots)\n'
id_clim_sst.close()



# Table to host 1 zonal profile per EXP:
vzc = nmp.zeros(nj0) ; # a zonal profile...


# Getting land-sea mask and coordinates:
bt.chck4f(vdic['MM_FILE'])
id_mask = Dataset(vdic['MM_FILE'])
xlon  = id_mask.variables['nav_lon'][:,:]
xlat  = id_mask.variables['nav_lat'][:,:]
imask = id_mask.variables['tmask'][0,:,:,:]
id_mask.close()






# Getting NEMO mean monthly climatology of temperature and salinity:
# ------------------------------------------------------------------

cvT3d=vdic['NN_T']
cvS3d=vdic['NN_S']
cf_nemo_mnmc = vdic['DIAG_D']+'/clim/mclim_'+CONFEXP+'_'+cy1+'-'+cy2+'_grid_T.nc4'

bt.chck4f(cf_nemo_mnmc)

id_nemo_mnmc = Dataset(cf_nemo_mnmc)

list_var = id_nemo_mnmc.variables.keys()
if vdic['NN_SST'] == cvT3d:
    SSTnemo = id_nemo_mnmc.variables[cvT3d][:,0,:,:]
else:
    SSTnemo = id_nemo_mnmc.variables[vdic['NN_SST']][:,:,:]

if vdic['NN_SSS'] == cvS3d:
    SSSnemo = id_nemo_mnmc.variables[cvS3d][:,0,:,:]
else:
    SSSnemo = id_nemo_mnmc.variables[vdic['NN_SSS']][:,:,:]

l_do_3d=True
if 'deptht' in list_var:
    vdepth = id_nemo_mnmc.variables['deptht'][:]
else:
    print 'WARNING: depth vector "deptht" not present in '+cf_nemo_mnmc+'!\n'
    l_do_3d=False

if cvT3d in list_var:
    Tnemo  = id_nemo_mnmc.variables[cvT3d][:,:,:,:]
    print '(has ',Tnemo.shape[0],' time snapshots)\n'
else:
    print 'WARNING: 3D NEMO T '+cvT3d+' not present in '+cf_nemo_mnmc+'!\n'
    l_do_3d=False

if cvS3d in list_var:
    Snemo  = id_nemo_mnmc.variables[cvS3d][:,:,:,:]
else:
    print 'WARNING: 3D NEMO S '+cvS3d+' not present in '+cf_nemo_mnmc+'!\n'
    l_do_3d=False

id_nemo_mnmc.close()

if l_do_3d:
    [ nt, nk, nj, ni ] = Tnemo.shape
    if (nk,nj,ni) != (nk0,nj0,ni0):
        print 'ERROR: 3D clim and NEMO file do no agree in shape!'
        print '       clim => '+str(ni0)+', '+str(nj0)+', '+str(nk0),' ('+vdic['F_T_CLIM_3D_12']+')'
        print '       NEMO => '+str(ni)+', '+str(nj)+', '+str(nk)
        sys.exit(0)
else:
    [ nt, nj, ni ] = SSTnemo.shape
    if (nj,ni) != (nj0,ni0):
        print 'ERROR: 3D clim and NEMO file do no agree in shape!'
        print '       clim => '+str(ni0)+', '+str(nj0),' ('+vdic['F_T_CLIM_3D_12']+')'
        print '       NEMO => '+str(ni)+', '+str(nj)
        sys.exit(0)




# Saving some array to avoid to call 'nmp.mean' all the time:

#Annual:
if l_do_3d:
    Tnemo_annual = nmp.zeros((nk,nj,ni))
    Tnemo_annual[:,:,:] = nmp.mean(Tnemo[:,:,:,:], axis=0)
    Snemo_annual = nmp.zeros((nk,nj,ni))
    Snemo_annual[:,:,:] = nmp.mean(Snemo[:,:,:,:], axis=0)
    Tclim_annual = nmp.zeros((nk,nj,ni))
    Tclim_annual[:,:,:] = nmp.mean(Tclim[:,:,:,:], axis=0)
Sclim_annual = nmp.zeros((nk0,nj,ni))
Sclim_annual[:,:,:] = nmp.mean(Sclim[:,:,:,:], axis=0)

SSTnemo_annual = nmp.zeros((nj,ni))
SSTnemo_annual[:,:] = nmp.mean(SSTnemo[:,:,:], axis=0)
SSSnemo_annual = nmp.zeros((nj,ni))
SSSnemo_annual[:,:] = nmp.mean(SSSnemo[:,:,:], axis=0)
SSTclim_annual = nmp.zeros((nj,ni))
SSTclim_annual[:,:] = nmp.mean(SSTclim[:,:,:], axis=0)

#JFM:
if l_do_3d:
    Tnemo_JFM = nmp.zeros((nk,nj,ni))
    Tnemo_JFM[:,:,:] = nmp.mean(Tnemo[:3,:,:,:], axis=0)
    Snemo_JFM = nmp.zeros((nk,nj,ni))
    Snemo_JFM[:,:,:] = nmp.mean(Snemo[:3,:,:,:], axis=0)
    Tclim_JFM = nmp.zeros((nk,nj,ni))
    Tclim_JFM[:,:,:] = nmp.mean(Tclim[:3,:,:,:], axis=0)
    Sclim_JFM = nmp.zeros((nk,nj,ni))
    Sclim_JFM[:,:,:] = nmp.mean(Sclim[:3,:,:,:], axis=0)

SSTnemo_JFM = nmp.zeros((nj,ni))
SSTnemo_JFM[:,:] = nmp.mean(SSTnemo[:3,:,:], axis=0)
SSSnemo_JFM = nmp.zeros((nj,ni))
SSSnemo_JFM[:,:] = nmp.mean(SSSnemo[:3,:,:], axis=0)
SSTclim_JFM = nmp.zeros((nj,ni))
SSTclim_JFM[:,:] = nmp.mean(SSTclim[:3,:,:], axis=0)

#JAS:
if l_do_3d:
    Tnemo_JAS = nmp.zeros((nk,nj,ni))
    Tnemo_JAS[:,:,:] = nmp.mean(Tnemo[6:9,:,:,:], axis=0)
    Snemo_JAS = nmp.zeros((nk,nj,ni))
    Snemo_JAS[:,:,:] = nmp.mean(Snemo[6:9,:,:,:], axis=0)
    Tclim_JAS = nmp.zeros((nk,nj,ni))
    Tclim_JAS[:,:,:] = nmp.mean(Tclim[6:9,:,:,:], axis=0)
    Sclim_JAS = nmp.zeros((nk,nj,ni))
    Sclim_JAS[:,:,:] = nmp.mean(Sclim[6:9,:,:,:], axis=0)

SSTnemo_JAS = nmp.zeros((nj,ni))
SSTnemo_JAS[:,:] = nmp.mean(SSTnemo[6:9,:,:], axis=0)
SSSnemo_JAS = nmp.zeros((nj,ni))
SSSnemo_JAS[:,:] = nmp.mean(SSSnemo[6:9,:,:], axis=0)
SSTclim_JAS = nmp.zeros((nj,ni))
SSTclim_JAS[:,:] = nmp.mean(SSTclim[6:9,:,:], axis=0)

if l_do_3d:
    jk100  = bt.find_index_from_value(100.  , vdepth) ; print 'jk100  = ', jk100,  '=> ', vdepth[jk100]
    jk1000 = bt.find_index_from_value(1000. , vdepth) ; print 'jk1000 = ', jk1000, '=> ', vdepth[jk1000]
    jk3000 = bt.find_index_from_value(3000. , vdepth) ; print 'jk3000 = ', jk3000, '=> ', vdepth[jk3000]
    tdj = [ jk100,   jk1000, jk3000  ]
    tdd_true = [ str(int(round(vdepth[jk100])))+'m' , str(int(round(vdepth[jk1000])))+'m' , str(int(round(vdepth[jk3000])))+'m' ]
    tdd      = [ '100m', '1000m', '3000m' ]
    print '\n', tdd_true[:], '\n'

# Creating 1D long. and lat.:
vlon = nmp.zeros(ni) ; vlon[:] = xlon[0,:]
ji_lat0 = nmp.argmax(xlat[nj-1,:])
vlat = nmp.zeros(nj) ; vlat[:] = xlat[:,ji_lat0]


# Time for figures:
# -----------------

if lfig0:

    if CC == 'CLIM':
        ctt = CONFEXP+': Mean Annual Zonal Anomaly of SST / Reynolds, ('+cy1+'-'+cy2+')'
    else:
        ctt = CONFEXP+': Mean Annual Zonal Anomaly of SST / '+CC+', ('+cy1+'-'+cy2+')'

    vzc[:] = bt.mk_zonal(SSTnemo_annual[:,:] - SSTclim_annual[:,:], imask[0,:,:])
    # Only at the end of all the experiments we do 2d plotting:
    bp.plot("zonal")(vlat, vzc, cfignm=path_fig+'1d_zonal_temp_anom_vs_'+CC, zmin=-5., zmax=5., dz=1.,
                     xmin=-75., xmax=65., czunit=r'$^{\circ}$C', cfig_type=fig_type,
                     ctitle=ctt)

    if CC == 'CLIM':
        ctt = CONFEXP+': Mean Annual Zonal Anomaly of SSS / WOA2009, ('+cy1+'-'+cy2+')'
    else:
        ctt = CONFEXP+': Mean Annual Zonal Anomaly of SSS / '+CC+', ('+cy1+'-'+cy2+')'

    vzc[:] = bt.mk_zonal(SSSnemo_annual[:,:] - Sclim_annual[0,:,:], imask[0,:,:])
    # Only at the end of all the experiments we do 2d plotting:
    bp.plot("zonal")(vlat, vzc, cfignm=path_fig+'1d_zonal_sali_anom_vs_'+CC , zmin=-2.5, zmax=2.5, dz=0.5,
                     xmin=-75., xmax=65., czunit='PSU', cfig_type=fig_type,
                     ctitle=ctt)



if lfig1:

    #                    SST / Reynolds
    # JFM
    bp.plot("2d")(vlon, vlat, SSTnemo_JFM[:,:] - SSTclim_JFM[:,:],
                  imask[0,:,:], tmin, tmax, dtemp,
                  corca=vdic['ORCA'], lkcont=False, cpal='RdBu_r', cfignm=path_fig+'dsst_JFM_'+CONFEXP+'_-_'+CC,
                  cbunit='K', cfig_type=fig_type,
                  ctitle='SST difference to '+CC+', JFM, '+CONFEXP+' ('+cy1+'-'+cy2+')',
                  lforce_lim=True)
    # JAS
    bp.plot("2d")(vlon, vlat, SSTnemo_JAS[:,:] - SSTclim_JAS[:,:],
                  imask[0,:,:], tmin, tmax, dtemp,
                  corca=vdic['ORCA'], lkcont=False, cpal='RdBu_r', cfignm=path_fig+'dsst_JAS_'+CONFEXP+'_-_'+CC,
                  cbunit='K', cfig_type=fig_type,
                  ctitle='SST difference to '+CC+', JAS, '+CONFEXP+' ('+cy1+'-'+cy2+')',
                  lforce_lim=True)

    # Annual
    bp.plot("2d")(vlon, vlat, SSTnemo_annual[:,:] - SSTclim_annual[:,:],
                  imask[0,:,:],  tmin, tmax, dtemp,
                  corca=vdic['ORCA'], lkcont=False, cpal='RdBu_r', cfignm=path_fig+'dsst_annual_'+CONFEXP+'_-_'+CC,
                  cbunit='K', cfig_type=fig_type,
                  ctitle='SST difference to '+CC+', '+CONFEXP+' ('+cy1+'-'+cy2+')',
                  lforce_lim=True)

    # Temperature 100m, 1000m... / climatology
    if l_do_3d:
        for jd in range(nmp.size(tdj)):
            jdepth = tdj[jd] ; cdepth = tdd[jd] ; cdepth_true = tdd_true[jd]    
            print '\n Treating depth '+str(vdepth[jdepth])+' !!!'
    
            if jd < 1:
                # JFM
                bp.plot("2d")(vlon, vlat, Tnemo_JFM[jdepth,:,:] - Tclim_JFM[jdepth,:,:],
                              imask[jdepth,:,:], tmin, tmax, dtemp,
                              corca=vdic['ORCA'], lkcont=False, cpal='RdBu_r', cfignm=path_fig+'dT_JFM_'+cdepth+'_'+CONFEXP+'_-_'+CC,
                              cbunit='K', cfig_type=fig_type,
                              ctitle='Temperature diff. to '+CC+' at '+cdepth_true+', JFM, '+CONFEXP+' ('+cy1+'-'+cy2+')',
                              lforce_lim=True)
                # JAS
                bp.plot("2d")(vlon, vlat, Tnemo_JAS[jdepth,:,:] - Tclim_JAS[jdepth,:,:],
                              imask[jdepth,:,:], tmin, tmax, dtemp,
                              corca=vdic['ORCA'], lkcont=False, cpal='RdBu_r', cfignm=path_fig+'dT_JAS_'+cdepth+'_'+CONFEXP+'_-_'+CC,
                              cbunit='K', cfig_type=fig_type,
                              ctitle='Temperature diff. to '+CC+' at '+cdepth_true+', JAS, '+CONFEXP+' ('+cy1+'-'+cy2+')',
                              lforce_lim=True)
    
            # Annual
            bp.plot("2d")(vlon, vlat, Tnemo_annual[jdepth,:,:] - Tclim_annual[jdepth,:,:],
                          imask[jdepth,:,:], tmin, tmax, dtemp,
                          corca=vdic['ORCA'], lkcont=False, cpal='RdBu_r', cfignm=path_fig+'dT_annual_'+cdepth+'_'+CONFEXP+'_-_'+CC,
                          cbunit='K', cfig_type=fig_type,
                          ctitle='Temperature diff. to '+CC+' at '+cdepth_true+', '+CONFEXP+' ('+cy1+'-'+cy2+')',
                          lforce_lim=True)






    #                   S S S
    # JFM
    bp.plot("2d")(vlon, vlat, SSSnemo_JFM[:,:] - Sclim_JFM[0,:,:],
                  imask[0,:,:], smin, smax, dsali, cpal='PiYG_r',
                  corca=vdic['ORCA'], lkcont=False, cfignm=path_fig+'dsss_JFM_'+CONFEXP+'_-_'+CC,
                  cbunit='PSU', cfig_type=fig_type,
                  ctitle='SSS difference to '+CC+', JFM, '+CONFEXP+' ('+cy1+'-'+cy2+')',
                  lforce_lim=True)
    # JAS
    bp.plot("2d")(vlon, vlat, SSSnemo_JAS[:,:] - Sclim_JAS[0,:,:],
                  imask[0,:,:], smin, smax, dsali, cpal='PiYG_r',
                  corca=vdic['ORCA'], lkcont=False, cfignm=path_fig+'dsss_JAS_'+CONFEXP+'_-_'+CC,
                  cbunit='PSU', cfig_type=fig_type,
                  ctitle='SSS difference to '+CC+', JAS, '+CONFEXP+' ('+cy1+'-'+cy2+')',
                  lforce_lim=True)
    # Annual
    bp.plot("2d")(vlon, vlat, SSSnemo_annual[:,:] - Sclim_annual[0,:,:],
                  imask[0,:,:], smin, smax, dsali, cpal='PiYG_r',
                  corca=vdic['ORCA'], lkcont=False, cfignm=path_fig+'dsss_annual_'+CONFEXP+'_-_'+CC,
                  cbunit='PSU', cfig_type=fig_type,
                  ctitle='SSS difference to '+CC+', '+CONFEXP+' ('+cy1+'-'+cy2+')',
                  lforce_lim=True)




    # Salinity 100m, 1000m... / climatology
    if l_do_3d:
        for jd in range(nmp.size(tdj)):
            jdepth = tdj[jd] ; cdepth = tdd[jd] ; cdepth_true = tdd_true[jd]    
            print '\n Treating depth '+str(vdepth[jdepth])+' !!!'
    
            if jd < 1:
                # JFM
                bp.plot("2d")(vlon, vlat, Snemo_JFM[jdepth,:,:] - Sclim_JFM[jdepth,:,:],
                              imask[jdepth,:,:], smin, smax, dsali, cpal='PiYG_r',
                              corca=vdic['ORCA'], lkcont=False, cfignm=path_fig+'dS_JFM_'+cdepth+'_'+CONFEXP+'_-_'+CC,
                              cbunit='PSU', cfig_type=fig_type,
                              ctitle='Salinity diff. to '+CC+' at '+cdepth_true+', JFM, '+CONFEXP+' ('+cy1+'-'+cy2+')',
                              lforce_lim=True)
                # JAS
                bp.plot("2d")(vlon, vlat, Snemo_JAS[jdepth,:,:] - Sclim_JAS[jdepth,:,:],
                              imask[jdepth,:,:], smin, smax, dsali, cpal='PiYG_r',
                              corca=vdic['ORCA'], lkcont=False, cfignm=path_fig+'dS_JAS_'+cdepth+'_'+CONFEXP+'_-_'+CC,
                              cbunit='PSU', cfig_type=fig_type,
                              ctitle='Salinity diff. to '+CC+' at '+cdepth_true+', JAS, '+CONFEXP+' ('+cy1+'-'+cy2+')',
                              lforce_lim=True)
    
            # Annual
            bp.plot("2d")(vlon, vlat, Snemo_annual[jdepth,:,:] - Sclim_annual[jdepth,:,:],
                          imask[jdepth,:,:], smin, smax, dsali, cpal='PiYG_r',
                          corca=vdic['ORCA'], lkcont=False, cfignm=path_fig+'dS_annual_'+cdepth+'_'+CONFEXP+'_-_'+CC,
                          cbunit='PSU', cfig_type=fig_type,
                          ctitle='Salinity diff. to '+CC+' at '+cdepth_true+', '+CONFEXP+' ('+cy1+'-'+cy2+')',
                          lforce_lim=True)
    
    
    
if lfig2 and i_do_sect==1 and l_do_3d : # Temperature and salinity for vertical sections

    vdico = bt.check_env_var(sys.argv[0], {'TS_SECTION_FILE'})

    vboxes, vlon1, vlat1, vlon2, vlat2, vTmin, vTmax, vSmin, vSmax = bt.read_coor(vdico['TS_SECTION_FILE'], ctype='float', lTS_bounds=True)
    
    js = 0
    for csname in vboxes:
        
        [ i1, i2, j1, j2 ] = bo.transect_zon_or_med(vlon1[js], vlon2[js], vlat1[js], vlat2[js], xlon, xlat)
        if i2>i1 and i2 < ni0-1: i2 = i2+1
        if j2>j1 and j2 < nj0-1: j2 = j2+1

        print '\n *** Section: '+csname+':'
        print ' lon1, lon2, lat1, lat2 =', vlon1[js], vlon2[js], vlat1[js], vlat2[js]
        print ' => i1, i2, j1, j2 =', i1, i2, j1, j2
        #print   ' xlon[j1,i1], xlon[j1,i2] =', xlon[j1,i1]-360., xlon[j1,i2]
        print ''

        if i1 > i2: print 'ERROR: temp_sal.py => i1 > i2 !'; sys.exit(0)
        if j1 > j2: print 'ERROR: temp_sal.py => j1 > j2 !'; sys.exit(0)

        #  N E M O
        #  ~~~~~~~
        #
        lzonal = False
        if i1 == i2:
            print ' ==> Meridional section!'
            vaxis = xlat[j1:j2+1,i1]
            ZT    = Tnemo_annual[:,j1:j2+1,i1] ; ZS = Snemo_annual[:,j1:j2+1,i1]
            OT    = Tclim_annual[:,j1:j2+1,i1] ; OS = Sclim_annual[:,j1:j2+1,i1]
            imsk  = imask[:,j1:j2+1,i1]
            cinfo = ', lon='+str(vlon1[js])
            xmn=vlat1[js]; xmx=vlat2[js]
        if j1 == j2:
            print ' ==> Zonal section!'            
            lzonal = True
            xmn=vlon1[js]; xmx=vlon2[js]
            vx = xlon[j1,i1:i2+1] ; vaxis = nmp.zeros(len(vx)) ; vaxis[:] = vx[:]
            if xmx<0. and xmn>0.:
                xmx = 360. + xmx                
            else:
                ivf = nmp.where(vx>180.); vaxis[ivf] = vx[ivf] - 360.
            ZT    = Tnemo_annual[:,j1,i1:i2+1] ; ZS = Snemo_annual[:,j1,i1:i2+1]
            OT    = Tclim_annual[:,j1,i1:i2+1] ; OS = Sclim_annual[:,j1,i1:i2+1]
            imsk  = imask[:,j1,i1:i2+1]
            cinfo = ', lat='+str(vlat1[js])            
            #if xmx < xmn and xmx <0.: xmx = 360. + xmx
            #print 'LOLO: xmn, xmx =>', xmn, xmx ; #sys.exit(0)

        Ta = vTmax[js] - vTmin[js]
        if Ta < 1. : print 'ERROR: temp_sal.py => Problem with your min and max for T for section "'+csname+'" !'; sys.exit(0)
        if Ta >= 25.:              dT = 1.
        if Ta >= 10. and Ta < 25.: dT = 0.5
        if Ta >=  5. and Ta < 10.: dT = 0.25
        if Ta >   0. and Ta <  5.: dT = 0.1

        Sa = vSmax[js] - vSmin[js]
        if Sa < 0.001 : print 'ERROR: temp_sal.py => Problem with your min and max for S for section "'+csname+'" !'; sys.exit(0)
        if Sa >= 3. :              dS = 0.1
        if Sa >= 1.5 and Sa < 3. : dS = 0.05
        if Sa >= 0.5 and Sa < 1.5: dS = 0.025
        if Sa >= 0.  and Sa < 0.5: dS = 0.01

        bp.plot("vert_section")(vaxis, vdepth, ZT, imsk, vTmin[js], vTmax[js], dT,
                                cpal='ncview_nrl', lzonal=lzonal, xmin=xmn, xmax=xmx, dx=5.,
                                cfignm=path_fig+'section_T_'+csname+'_'+CONFEXP, cbunit=r'$^{\circ}$C', cxunit=r'Latitude ($^{\circ}$N)',
                                czunit='Depth (m)', ctitle='Temperature, ('+cy1+'-'+cy2+'), '+csname+', '+CONFEXP+cinfo,
                                cfig_type=fig_type, lforce_lim=False, i_cb_subsamp=2)
        
        bp.plot("vert_section")(vaxis, vdepth, ZS, imsk, vSmin[js], vSmax[js], dS,
                                cpal='ncview_jaisnb', lzonal=lzonal, xmin=xmn, xmax=xmx, dx=5.,
                                cfignm=path_fig+'section_S_'+csname+'_'+CONFEXP, cbunit='PSU', cxunit=r'Latitude ($^{\circ}$N)',
                                czunit='Depth (m)', ctitle='Salinity, ('+cy1+'-'+cy2+'), '+csname+', '+CONFEXP+cinfo,
                                cfig_type=fig_type, lforce_lim=False, i_cb_subsamp=2)

        #  OBS:
        bp.plot("vert_section")(vaxis, vdepth, OT, imsk, vTmin[js], vTmax[js], dT,
                                cpal='ncview_nrl', lzonal=lzonal, xmin=xmn, xmax=xmx, dx=5.,
                                cfignm=path_fig+'section_T_'+csname+'_'+CC, cbunit=r'$^{\circ}$C',
                                cxunit=r'Latitude ($^{\circ}$N)',
                                czunit='Depth (m)', ctitle='Temperature, '+csname+', '+CC+cinfo,
                                cfig_type=fig_type, lforce_lim=False, i_cb_subsamp=2)
        #
        bp.plot("vert_section")(vaxis, vdepth, OS, imsk, vSmin[js], vSmax[js], dS,
                                cpal='ncview_jaisnb', lzonal=lzonal, xmin=xmn, xmax=xmx, dx=5.,
                                cfignm=path_fig+'section_S_'+csname+'_'+CC, cbunit='PSU',
                                cxunit=r'Latitude ($^{\circ}$N)',
                                czunit='Depth (m)', ctitle='Salinity, '+csname+', '+CC+cinfo,
                                cfig_type=fig_type, lforce_lim=False, i_cb_subsamp=2)

        js=js+1


print '\n Bye!'
