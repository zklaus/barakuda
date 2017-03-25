#!/usr/bin/env python
#
#       B a r a K u d a
#
#    Spatially integrates surface fluxes (on NEMO grid) on all the 2D regions
#    defined in the basin_mask file
#
#       L. Brodeau, March 2017
#

import sys
import os
import numpy as nmp

from netCDF4 import Dataset
from string  import replace

import barakuda_tool as bt
import barakuda_orca as bo
import barakuda_ncio as bnc

venv_needed = {'ORCA','EXP','DIAG_D','MM_FILE','BM_FILE','FILE_FLX_SUFFIX','TSTAMP',
               'NEMO_SAVED_FILES','NN_RNF','NN_FWF','NN_P','NN_CLV','NN_E'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

cdir_out = vdic['DIAG_D']+'/flux_int_basins'
os.system('mkdir -p '+cdir_out)

cFgrid = vdic['FILE_FLX_SUFFIX'] ; # suffix of files containing surface fluxes!

CONFEXP = vdic['ORCA']+'-'+vdic['EXP']

if len(sys.argv) != 3:
    print 'Usage : sys.argv[1] <ORCA1_EXP_'+cFgrid+'.nc> <year>'; sys.exit(0)

cnexec  = sys.argv[0]
cf_flx_in = sys.argv[1]
cyear   = sys.argv[2] ; jyear = int(cyear); cyear = '%4.4i'%jyear

print 'Current year is '+cyear+' !\n'


# Checking if the land-sea mask file is here:
for cf in [vdic['MM_FILE'], vdic['BM_FILE']]:
    if not os.path.exists(cf):
        print 'Mask file '+cf+' not found'; sys.exit(0)

# Reading the grid metrics:
id_mm = Dataset(vdic['MM_FILE'])
list_variables = id_mm.variables.keys()
mask_glo = id_mm.variables['tmask'][0,0,:,:]
xe1t     = id_mm.variables['e1t'][0,:,:]
xe2t     = id_mm.variables['e2t'][0,:,:]
id_mm.close()

(nj, ni) = mask_glo.shape

# Grid cell area:
Xa   = nmp.zeros((nj, ni))
Xa[:,:] = 1.E-6*xe1t[:,:]*xe2t[:,:]*mask_glo[:,:]
del xe1t, xe2t


print ' *** Global ocean surface =', 1.E-6*nmp.sum(Xa[:-2,:-2]), '10^6 km^2\n'


print 'Opening different basin masks in file '+vdic['BM_FILE']
list_basin_names, list_basin_lgnms = bo.get_basin_info(vdic['BM_FILE'])
nb_basins = len(list_basin_names) ; # Global is #0, it is included!
mask      = nmp.zeros((nb_basins,nj,ni))
msk_tmp   = nmp.zeros((nj,ni))


# Geting basin masks and their zooming indices:
i1 = nmp.zeros(nb_basins, dtype=nmp.int)
i2 = nmp.zeros(nb_basins, dtype=nmp.int)
j1 = nmp.zeros(nb_basins, dtype=nmp.int)
j2 = nmp.zeros(nb_basins, dtype=nmp.int)

mask[0,:,:] = mask_glo[:,:] ; # global
id_bm = Dataset(vdic['BM_FILE'])
for jb in range(nb_basins) :
    if jb > 0:
        msk_tmp[:,:] = id_bm.variables['tmask'+list_basin_names[jb]][:,:]
        mask[jb,:,:] = msk_tmp[:,:]*mask_glo[:,:]    
    # Decrasing the domain size if possible:
    (i1[jb],j1[jb],i2[jb],j2[jb]) = bo.shrink_domain(mask[jb,:,:])
    #(vjj , vji)  = nmp.where(mask[jb,:,:]>0.5)
    #j1[jb] = max( nmp.min(vjj)-2 , 0    )
    #i1[jb] = max( nmp.min(vji)-2 , 0    )
    #j2[jb] = min( nmp.max(vjj)+2 , nj-1 ) + 1
    #i2[jb] = min( nmp.max(vji)+2 , ni-1 ) + 1            
    #if (i1[jb],i2[jb]) == (0,ni): i2[jb] = i2[jb]-2 ; # Mind east-west periodicity overlap of 2 points...
id_bm.close()

del mask_glo, msk_tmp


# Getting list of variables available in input NEMO file:
id_in = Dataset(cf_flx_in)
list_variables = id_in.variables.keys()
id_in.close()




# Fresh water fluxes (expected in kg/m2/s == mm/s!):
rmult = 1.E-3 ; # to Sverdrup!
cunit='Sv'

vvar = [ vdic['NN_FWF'], vdic['NN_E'], vdic['NN_P'], vdic['NN_RNF'], vdic['NN_CLV'] ]
vnnm = [ 'EmPmR'       ,   'E'       ,        'P'  ,       'R'     ,   'ICalv'      ]
nbv  = len(vvar)

# Array contening results for freshwater fluxes:
v_fwf_m = nmp.zeros((12,nb_basins,nbv))

jv = 0
for cvar in vvar:

    if cvar == 'X' or (not cvar in list_variables):
        print '\n WARNING ('+cnexec+'): skipping '+vnnm[jv]+' = "'+cvar+'"!'
        if cvar != 'X': print '     => because not found in '+cf_flx_in+' !\n'
    else:
        print '\n +++ '+cnexec+' => treating '+vnnm[jv]+' ('+cvar+') !'
        print '      ==> variable '+cvar
        id_in = Dataset(cf_flx_in)
        Xd_m = id_in.variables[cvar][:,:,:]
        cu = id_in.variables[cvar].units
        id_in.close()
        (Nt,nj0,ni0) = nmp.shape(Xd_m)
        print '      ==> variable '+cvar+' read ! ('+cu+') =>', (Nt,nj0,ni0)
        
        if Nt != 12: print 'ERROR: '+cnexec+' => only treating monthly data so far...'; sys.exit(0)
    
        if (nj0,ni0) != (nj,ni): print 'ERROR: '+cnexec+' => Field and metrics do not agree in size!'; sys.exit(0)
    
        Xd_m_x_A = nmp.zeros((Nt,nj,ni))
        for jt in range(Nt):
            Xd_m_x_A[jt,:,:] = Xd_m[jt,:,:]*Xa[:,:]
            
        jb = 0
        for cbasin in list_basin_names[:]:
            colnm = list_basin_lgnms[jb]
            print '   ===> '+cvar+' in basin '+cbasin+' ('+colnm+')'
            if (i1[jb],j1[jb],i2[jb],j2[jb]) != (0,0,ni,nj): print '       ===> zooming on i1,j1 -> i2,j2:', i1[jb],j1[jb],'->',i2[jb],j2[jb]    
    
            for jt in range(Nt):
                v_fwf_m[jt,jb,jv] = nmp.sum( rmult*Xd_m_x_A[jt,j1[jb]:j2[jb],i1[jb]:i2[jb]] * mask[jb,j1[jb]:j2[jb],i1[jb]:i2[jb]] )
    
            jb = jb + 1 ; # next basin!

        print '\n +++ '+cnexec+' => Done with 2D-summing of variable '+cvar+'!\n'
        del Xd_m
        
    jv = jv + 1 ; # next variable!


# Time to write netcdf file:

vtime = nmp.zeros(12)
for jt in range(12): vtime[jt] = float(jyear) + (float(jt)+0.5)*1./float(12)
#print ' * Calendar: ', vtime[:]


c1 = '2D-integral of '
jb = 0
for cbasin in list_basin_names[:]:
    c2 = ') on region '+list_basin_lgnms[jb]
    cf_out = cdir_out+'/fwf_int_'+CONFEXP+'_'+cbasin+'.nc'
    bnc.wrt_appnd_1d_series(vtime, v_fwf_m[:,jb,0], cf_out, vnnm[0],
                            cu_t='year', cu_d=cunit,            cln_d =c1+vnnm[0]+' ('+vvar[0]+c2,
                            vd2=v_fwf_m[:,jb,1], cvar2=vnnm[1], cln_d2=c1+vnnm[1]+' ('+vvar[1]+c2,
                            vd3=v_fwf_m[:,jb,2], cvar3=vnnm[2], cln_d3=c1+vnnm[2]+' ('+vvar[2]+c2,
                            vd4=v_fwf_m[:,jb,3], cvar4=vnnm[3], cln_d4=c1+vnnm[3]+' ('+vvar[3]+c2,
                            vd5=v_fwf_m[:,jb,4], cvar5=vnnm[4], cln_d5=c1+vnnm[4]+' ('+vvar[4]+c2)
    jb = jb + 1


print '\n *** EXITING '+cnexec+' for year '+cyear+' and variable '+cvar+'!\n'
