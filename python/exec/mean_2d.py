#!/usr/bin/env python
#
#       B a r a K u d a
#
#     Generate misc. spatial 2D averaging out of NEMO output files...
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



# Box nino 3.4:
lon1_nino = 360. - 170.  ; # east
lat1_nino = -5.
lon2_nino = 360. - 120.  ; # east
lat2_nino = 5.


#cv_evb = 'evap_ao_cea' ; # debug evap in ec-earth...

venv_needed = {'ORCA','EXP','DIAG_D','MM_FILE','BM_FILE','NEMO_SAVED_FILES','FILE_FLX_SUFFIX',
               'NN_SST','NN_SSS','NN_SSH','NN_MLD','TSTAMP','NN_FWF','NN_EMP','NN_P','NN_RNF',
               'NN_CLV','NN_E','NN_QNET','NN_QSOL'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFEXP = vdic['ORCA']+'-'+vdic['EXP']

if len(sys.argv) != 3:
    print 'Usage : sys.argv[1] <ORCA1_EXP_grid_T.nc> <year>'
    sys.exit(0)

cnexec  = sys.argv[0]
cf_T_in = sys.argv[1]
cyear   = sys.argv[2] ; jyear = int(cyear); cyear = '%4.4i'%jyear

print 'Current year is '+cyear+' !\n'

vtime = nmp.zeros(12)
for jt in range(12):
    vtime[jt] = float(jyear) + (float(jt)+0.5)*1./12.




# Checking if the land-sea mask file is here:
for cf in [vdic['MM_FILE'], vdic['BM_FILE']]:
    if not os.path.exists(cf):
        print 'Mask file '+cf+' not found'; sys.exit(0)

# Reading the grid metrics:
id_mm = Dataset(vdic['MM_FILE'])
list_variables = id_mm.variables.keys()
rmask  = id_mm.variables['tmask'][0,0,:,:]
xlon   = id_mm.variables['glamt'][0,:,:]
xlat   = id_mm.variables['gphit'][0,:,:]
xe1t   = id_mm.variables['e1t'][0,:,:]
xe2t   = id_mm.variables['e2t'][0,:,:]
id_mm.close()

[ nj, ni ] = rmask.shape


# About heat and freshwater fluxes:
cfe_sflx = vdic['FILE_FLX_SUFFIX']
l_fwf = False
l_htf = False
if cfe_sflx in vdic['NEMO_SAVED_FILES']:
    if vdic['NN_FWF']  != 'X': l_fwf = True
    if vdic['NN_QNET'] != 'X': l_htf = True
    cf_F_in = replace(cf_T_in, 'grid_T', cfe_sflx)


Xa   = nmp.zeros((nj, ni))
Xa[:,:] = xe1t[:,:]*xe2t[:,:]
del xe1t, xe2t

Xarea_t = nmp.zeros((nj, ni))
Xarea_t[:,:] = Xa[:,:]*rmask[:,:]*1.E-6 ; # [10^6 m^2] !
Socean = nmp.sum( Xarea_t[:-2,:-2] ) ; # Mind the 2 point folding overlap, must not be included in sum!
print '\n  *** Surface of the ocean = ', Socean* 1.E-6, '  [10^6 km^2]\n'



print 'Reading all basin masks in file '+vdic['BM_FILE']

list_basin_names, list_basin_lgnms = bo.get_basin_info(vdic['BM_FILE'])
nb_basins = len(list_basin_names)
mask = nmp.zeros((nb_basins,nj,ni))
msk_tmp = nmp.zeros((nj,ni))
mask[0,:,:] = rmask[:,:] ; # global

id_bm = Dataset(vdic['BM_FILE'])
for jb in range(1,nb_basins) :
    msk_tmp[:,:] = id_bm.variables['tmask'+list_basin_names[jb]][:,:]
    mask[jb,:,:] = msk_tmp[:,:]*rmask[:,:]
id_bm.close()

del rmask, msk_tmp





#######################################################################
# Time-series of globally averaged surface heat and freshwater fluxes #
######################################################################

# Heat fluxes
if l_htf:

    print '\n\n +++ '+cnexec+' => Starting heat flux diags!'

    cv_qnt = vdic['NN_QNET']
    cv_qsr = vdic['NN_QSOL']

    id_in = Dataset(cf_F_in)
    list_variables = id_in.variables.keys()
    l_qnt = False
    if  cv_qnt in list_variables[:]:
        l_qnt = True
        QNT_m = id_in.variables[cv_qnt][:,:,:]
        print '   *** Qnet ('+cv_qnt+') read!'             
    l_qsr = False
    if  cv_qsr in list_variables[:]:
        l_qsr = True
        QSR_m = id_in.variables[cv_qsr][:,:,:]
        print '   *** Qsol ('+cv_qsr+') read!'
    id_in.close()
    l_htf = l_qnt ; # Forgeting heat flux if both Qnet is missing...

    [ nt, nj0, ni0 ] = QNT_m.shape
    vtime = nmp.zeros(nt)
    for jt in range(nt): vtime[jt] = float(jyear) + (float(jt)+0.5)*1./12.
    vqnt = [] ; vqsr = []
    if l_qnt: vqnt = nmp.zeros(nt)
    if l_qsr: vqsr = nmp.zeros(nt)
    for jt in range(nt):
        # Mind the 2 point folding overlap, must not be included in sum!
        if l_qnt: vqnt[jt] = nmp.sum( QNT_m[jt,:-2,:-2]*Xarea_t[:-2,:-2] ) * 1.E-9 ;  # to PW
        if l_qsr: vqsr[jt] = nmp.sum( QSR_m[jt,:-2,:-2]*Xarea_t[:-2,:-2] ) * 1.E-9 ;  # to PW

    cf_out   = vdic['DIAG_D']+'/mean_htf_'+CONFEXP+'_GLO.nc'
    bnc.wrt_appnd_1d_series(vtime, vqnt, cf_out, 'Qnet',
                            cu_t='year', cu_d='PW', cln_d ='Globally averaged net heat flux (nemo:'+cv_qnt+')',
                            vd2=vqsr, cvar2='Qsol', cln_d2='Globally averaged net solar heat flux (nemo:'+cv_qsr+')')
    print ' +++ '+cnexec+' => Done with heat flux diags!\n'

# Freshwater fluxes
if l_fwf:

    print '\n\n +++ '+cnexec+' => Starting freshwater flux diags!'

    cv_fwf = vdic['NN_FWF']
    cv_emp = vdic['NN_EMP']
    cv_prc = vdic['NN_P']
    cv_rnf = vdic['NN_RNF']
    cv_clv = vdic['NN_CLV']
    cv_evp = vdic['NN_E']
    # cv_evb (top of file...)

    
    id_in = Dataset(cf_F_in)
    list_variables = id_in.variables.keys()

    if not cv_fwf in list_variables[:]: print 'PROBLEM with fwf! '+cnexec+'!'; sys.exit(0)

    FWF_m = id_in.variables[cv_fwf][:,:,:]
    print '   *** E-P-R ('+cv_fwf+') read!'

    l_emp = False
    if  cv_emp in list_variables[:]:
        l_emp = True
        EMP_m = id_in.variables[cv_emp][:,:,:]
        print '   *** E-P ('+cv_emp+') read!'
             
    l_prc = False
    if  cv_prc in list_variables[:]:
        l_prc = True
        PRC_m = id_in.variables[cv_prc][:,:,:]
        print '   *** P ('+cv_prc+') read!'

    l_rnf = False
    if  cv_rnf in list_variables[:]:
        l_rnf = True
        RNF_m = id_in.variables[cv_rnf][:,:,:]
        print '   *** Runoffs ('+cv_rnf+') read!'

    l_clv = False
    if  cv_clv in list_variables[:]:
        l_clv = True
        CLV_m = id_in.variables[cv_clv][:,:,:]
        print '   *** Calving ('+cv_clv+') read!'

    l_evp = False
    if  cv_evp in list_variables[:]:
        l_evp = True
        EVP_m = id_in.variables[cv_evp][:,:,:]
        print '   *** Calving ('+cv_evp+') read!'

    #l_evb = False
    #if  cv_evb in list_variables[:]:
    #    l_evb = True
    #    EVB_m = id_in.variables[cv_evb][:,:,:]
    #    print '   *** Calving ('+cv_evb+') read!'

    id_in.close()

               
    [ nt, nj0, ni0 ] = FWF_m.shape
    
    if l_emp and not l_rnf:
        l_rnf = True
        RNF_m = nmp.zeros((nj0,ni0))
        RNF_m = - ( FWF_m - EMP_m )

    vtime = nmp.zeros(nt)
    for jt in range(nt): vtime[jt] = float(jyear) + (float(jt)+0.5)*1./12.

    vfwf = nmp.zeros(nt)
    
    vemp = [] ; vrnf = [] ; vprc = [] ; vclv = [] ; vevp = [] ; #vevb = []
    if l_emp: vemp = nmp.zeros(nt)
    if l_rnf: vrnf = nmp.zeros(nt)
    if l_prc: vprc = nmp.zeros(nt)
    if l_clv: vclv = nmp.zeros(nt)
    if l_evp: vevp = nmp.zeros(nt)
    #if l_evb: vevb = nmp.zeros(nt)


    for jt in range(nt):
        # Mind the 2 point folding overlap, must not be included in sum!
        vfwf[jt]           = nmp.sum( FWF_m[jt,:-2,:-2]*Xarea_t[:-2,:-2] ) * 1.E-3 ;  # to Sv
        if l_emp: vemp[jt] = nmp.sum( EMP_m[jt,:-2,:-2]*Xarea_t[:-2,:-2] ) * 1.E-3 ;  # to Sv
        if l_rnf: vrnf[jt] = nmp.sum( RNF_m[jt,:-2,:-2]*Xarea_t[:-2,:-2] ) * 1.E-3 ;  # to Sv
        if l_prc: vprc[jt] = nmp.sum( PRC_m[jt,:-2,:-2]*Xarea_t[:-2,:-2] ) * 1.E-3 ;  # to Sv
        if l_clv: vclv[jt] = nmp.sum( CLV_m[jt,:-2,:-2]*Xarea_t[:-2,:-2] ) * 1.E-3 ;  # to Sv
        if l_evp: vevp[jt] = nmp.sum( EVP_m[jt,:-2,:-2]*Xarea_t[:-2,:-2] ) * 1.E-3 ;  # to Sv
        #if l_evb: vevb[jt] = nmp.sum( EVB_m[jt,:-2,:-2]*Xarea_t[:-2,:-2] ) * 1.E-3 ;  # to Sv

    cf_out   = vdic['DIAG_D']+'/mean_fwf_'+CONFEXP+'_GLO.nc'

    bnc.wrt_appnd_1d_series(vtime, vfwf, cf_out, 'EmPmR',
                            cu_t='year', cu_d='Sv',  cln_d ='Globally averaged net freshwater flux (nemo:'+cv_fwf+')',
                            vd2=vemp, cvar2='EmP',   cln_d2='Globally averaged Evap - Precip (nemo:'+cv_emp+')',
                            vd3=vrnf, cvar3='R',     cln_d3='Globally averaged continental runoffs',
                            vd4=vprc, cvar4='P',     cln_d4='Globally averaged total precip (nemo:'+cv_prc+')',
                            vd5=vclv, cvar5='ICalv', cln_d5='Globally averaged ice calving from icebergs (nemo:'+cv_clv+')',
                            vd6=vevp, cvar6='E',     cln_d6='Globally averaged evaporation (nemo:'+cv_evp+')')
    #
    #                            vd7=vevb, cvar7='Eb',    cln_d7='Globally averaged evaporation with sea-ice consideration (nemo:'+cv_evb+')'
    #                        )

    print ' +++ '+cnexec+' => Done with freshwater flux diags!\n'






####################################
# MLD time serie in different boxes:
####################################

#print '\n\n +++ '+cnexec+' => Starting MLD diags!'
#l_mld = False
#print '\nSpatially-averaged MLD in different boxes'
#
#cvar = vdic['NN_MLD']
#id_in = Dataset(cf_T_in)
#list_variables = id_in.variables.keys()
#if cvar in list_variables[:]: # check if MLD variable is present!
#    MLD_m = id_in.variables[cvar][:,:,:]
#    print '   *** MLD ('+cvar+') found and read!'
#    l_mld = True
#else:
#    print '   *** OOPS! MLD ('+cvar+') not found, skipping MLD time series diag...'
#id_in.close()

#if l_mld:
#
#    [ nt, nj0, ni0 ] = MLD_m.shape
#
#    if nt != 12: print 'ERROR: '+cnexec+' => only treating monthly data so far...'; sys.exit(0)
#
#    if [ nj0, ni0 ] != [ nj, ni ]: print 'ERROR: '+cnexec+' => Field and metrics do not agree in size!'; sys.exit(0)
#
#    vtime = nmp.zeros(nt)
#    for jt in range(nt): vtime[jt] = float(jyear) + (float(jt)+0.5)*1./12.
#
#    mask2d = nmp.zeros((nj,ni))
#
#    # Reading boxes definitions into barakuda_orca.py:
#    cname_b = bo.cname_mld_boxes
#    nb_boxes = len(cname_b)
#
#    for ib in range(nb_boxes):
#
#        cbox = cname_b[ib] ; print '    *** treating '+cvar+' for '+cbox+', ('+bo.clgnm_mld_boxes[ib]+')'
#
#        i1 = nmp.argmax(xlat[nj-1,:])
#        j1 = 0 ; i2 = ni-1 ; j2 = nj-1
#
#        rx1 = bo.r_lon_p1_mld[ib] ; rx2 = bo.r_lon_p2_mld[ib] ; ry1 = bo.r_lat_p1_mld[ib] ; ry2 = bo.r_lat_p2_mld[ib]
#
#        # Need to itterate because ORCA grid distorded in the North...
#        vold = [ -999, -999, -999, -999 ] ;  itt = 0
#        while [ i1, i2, j1, j2 ] != vold and itt < 10 :
#            itt = itt+1
#            #print ' .... itt =', itt
#            vold = [ i1, i2, j1, j2 ]
#            #print 'seraching for rx1, rx2, ry1, ry2 = ', rx1, rx2, ry1, ry2
#            if ry1 > -900.: j1 = bt.find_index_from_value( ry1, xlat[:,i1] )
#            if rx1 > -900.: i1 = bt.find_index_from_value( rx1, xlon[j1,:] )
#            if rx2 > -900.: i2 = bt.find_index_from_value( rx2, xlon[j2,:] )
#            if ry1 > -900.: j1 = bt.find_index_from_value( ry1, xlat[:,i1] )
#            if ry2 > -900.: j2 = bt.find_index_from_value( ry2, xlat[:,i2] )
#            #print '   => i1, i2, j1, j2 =>', i1, i2, j1, j2, '\n'
#
#
#
#        mask2d[:,:] = 0.
#        mask2d[j1:j2,i1:i2] = mask[0,j1:j2,i1:i2]
#
#        Vts = bo.mean_2d(MLD_m, mask2d[:,:], Xa[:,:])
#
#        # NETCDF:
#        cf_out   = vdic['DIAG_D']+'/mean_'+cvar+'_'+CONFEXP+'_'+cbox+'.nc' ;  cv1 = cvar
#
#        bnc.wrt_appnd_1d_series(vtime, Vts, cf_out, cv1,
#                                cu_t='year', cu_d='m', cln_d='2D-average of '+cvar+' on rectangular box '+cbox)
#
#print ' +++ '+cnexec+' => Done with MLD diags!\n'


###############################################
# 2D (surface) averaging for T,S, SSH and MLD #
###############################################

print '\n\n +++ '+cnexec+' => Starting 2D surface averaging diags!'

id_in = Dataset(cf_T_in)
list_variables = id_in.variables.keys()
jvar = 0
for cvar in [ vdic['NN_SST'], vdic['NN_SSS'], vdic['NN_SSH'], vdic['NN_MLD'] ]:
    
    if cvar in list_variables:
        
        print '  *** reading '+cvar+' into '+cf_T_in
        if cvar in ['thetao','so','votemper','vosaline']:
            Xs_m = id_in.variables[cvar][:,0,:,:]
        else:
            Xs_m = id_in.variables[cvar][:,:,:]
        print '  ...read!'

        ( nt, nj0, ni0 ) = Xs_m.shape

        if nt != 12: print 'ERROR: '+cnexec+' => only treating monthly data so far...'; sys.exit(0)

        if ( nj0, ni0 ) != ( nj, ni ): print 'ERROR: '+cnexec+' => Field and metrics do not agree in size!'; sys.exit(0)

        if jvar == 0:
            vtime = nmp.zeros(nt)
            for jt in range(nt): vtime[jt] = float(jyear) + (float(jt)+0.5)*1./12.
            print ' * Montly calendar: ', vtime[:]

        joce = 0
        for cocean in list_basin_names[:]:

            print 'Treating '+cvar+' for '+cocean

            Vts = bo.mean_2d(Xs_m, mask[joce,:,:], Xa[:,:])

            cf_out   = vdic['DIAG_D']+'/mean_'+cvar+'_'+CONFEXP+'_'+cocean+'.nc' ;  cv1 = cvar
            bnc.wrt_appnd_1d_series(vtime, Vts, cf_out, cv1,
                                    cu_t='year', cu_d='m', cln_d='2D-average of '+cvar+' over '+list_basin_lgnms[joce])

            joce = joce + 1
        jvar = jvar + 1
        
    else:
        print 'WARNING: '+cnexec+' => '+cvar+' was not found in '+cf_T_in

id_in.close()
print ' +++ '+cnexec+' => Done with 2D surface averaging diags!\n'






###################
# El nino box 3.4 #
###################

print '\n\n +++ '+cnexec+' => Starting ENSO diag!'

print lon1_nino, lon2_nino, lat1_nino, lat2_nino
( i1, j1 ) = bo.ij_from_xy(lon1_nino, lat1_nino, xlon, xlat)
( i2, j2 ) = bo.ij_from_xy(lon2_nino, lat2_nino, xlon, xlat)
print ' Nino box 3.4, longitude: '+str(xlon[j1,i1])+' => '+str(xlon[j2,i2])+' \ latitude: '+str(xlat[j1,i1])+' => '+str(xlat[j2,i2])

id_in = Dataset(cf_T_in)
if vdic['NN_SST'] == 'thetao':
    Xs_m = id_in.variables[vdic['NN_SST']][:,0,:,:]
else:
    Xs_m = id_in.variables[vdic['NN_SST']][:,:,:]
id_in.close()

Vts = bo.mean_2d(Xs_m[:,j1:j2+1,i1:i2+1], mask[0,j1:j2+1,i1:i2+1], Xa[j1:j2+1,i1:i2+1])

cunit='deg.C'
if Vts[0] > 100.: cunit='K'

cf_out   = vdic['DIAG_D']+'/Nino34_'+CONFEXP+'.nc' ;  cv1 = vdic['NN_SST']
bnc.wrt_appnd_1d_series(vtime, Vts, cf_out, cv1, cu_t='year', cu_d=cunit,
                        cln_d='2D-average of SST over Nino box 3.4 for ENSO computation later')

print ' +++ '+cnexec+' => Done with ENSO diags!\n'






##########################################################
# AMO (Atlantic Multidecadal Oscillation) Index
##########################################################

print '\n\n +++ '+cnexec+' => Starting AMO diag!'
# Finding index jb of Atlantic mask:
jb=0
while list_basin_names[jb] != 'atl' : jb=jb+1    
print ' *** index Atl is ', jb, list_basin_names[jb]
msk_natl = nmp.zeros((nj,ni))
msk_natl[:,:] = mask[jb,:,:]
msk_natl[nmp.where(xlat< 2.)] = 0
msk_natl[nmp.where(xlat>68.)] = 0
#debug:bnc.write_2d_mask('mask_AMO.nc', msk_natl)

id_in = Dataset(cf_T_in)
if vdic['NN_SST'] == 'thetao':
    Xs_m = id_in.variables[vdic['NN_SST']][:,0,:,:]
else:
    Xs_m = id_in.variables[vdic['NN_SST']][:,:,:]
id_in.close()

Vts = bo.mean_2d(Xs_m, msk_natl, Xa)

cunit='deg.C'
if Vts[0] > 100.: cunit='K'

cf_out   = vdic['DIAG_D']+'/mean_SST_NAtl_'+CONFEXP+'.nc' ;  cv1 = vdic['NN_SST']
bnc.wrt_appnd_1d_series(vtime, Vts, cf_out, cv1, cu_t='year', cu_d=cunit,
                        cln_d='2D-average of SST over North Atlantic to compute AMO later')
del msk_natl, Vts
print ' +++ '+cnexec+' => Done with AMO diag!\n'










print '\n *** EXITING '+cnexec+' for year '+cyear+' !\n'
