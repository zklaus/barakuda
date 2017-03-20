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
rmask  = id_mm.variables['tmask'][0,0,:,:]
xe1t   = id_mm.variables['e1t'][0,:,:]
xe2t   = id_mm.variables['e2t'][0,:,:]
id_mm.close()

(nj, ni) = rmask.shape

# Grid cell area:
Xa   = nmp.zeros((nj, ni))
Xa[:,:] = 1.E-6*xe1t[:,:]*xe2t[:,:]*rmask[:,:]
del xe1t, xe2t


print ' *** Global ocean surface =', 1.E-6*nmp.sum(Xa[:-2,:-2]), '10^6 km^2\n'


print 'Opening different basin masks in file '+vdic['BM_FILE']
list_basin_names, list_basin_lgnms = bo.get_basin_info(vdic['BM_FILE'])
nb_basins = len(list_basin_names)
mask      = nmp.zeros((nb_basins,nj,ni))
msk_tmp   = nmp.zeros((nj,ni))
mask[0,:,:] = rmask[:,:] ; # global

id_bm = Dataset(vdic['BM_FILE'])
for joce in range(1,nb_basins) :
    msk_tmp[:,:] = id_bm.variables['tmask'+list_basin_names[joce]][:,:]
    mask[joce,:,:] = msk_tmp[:,:]*rmask[:,:]
id_bm.close()

del rmask, msk_tmp


# Getting list of variables available in input NEMO file:
id_in = Dataset(cf_flx_in)
list_variables = id_in.variables.keys()
id_in.close()


# Fresh water fluxes (expected in kg/m2/s == mm/s!):
rmult = 1.E-3 ; # to Sverdrup!
cunit='Sv'

# 1 RUNOFF !
vvar = [ vdic['NN_FWF'], vdic['NN_P'], vdic['NN_RNF'], vdic['NN_CLV'], vdic['NN_E'] ]

for cvar in vvar:

    if cvar == 'X' or (not cvar in list_variables):
        print '\n WARNING ('+cnexec+'): skipping NN_RNF="'+cvar+'"!'
        if cvar != 'X': print '     => because not found in '+cf_flx_in+' !\n'
    else:
        print '\n +++ '+cnexec+' => treating runoff!'
        print '      ==> variable '+cvar
        id_in = Dataset(cf_flx_in)
        Xd_m = id_in.variables[cvar][:,:,:]
        cu = id_in.variables[cvar].units
        id_in.close()
        (nt,nj0,ni0) = nmp.shape(Xd_m)
        print '      ==> variable '+cvar+' read ! ('+cu+') =>', (nt,nj0,ni0)
        
        if nt != 12: print 'ERROR: '+cnexec+' => only treating monthly data so far...'; sys.exit(0)
    
        if (nj0,ni0) != (nj,ni): print 'ERROR: '+cnexec+' => Field and metrics do not agree in size!'; sys.exit(0)
    
        vtime = nmp.zeros(nt)
        for jt in range(nt): vtime[jt] = float(jyear) + (float(jt)+0.5)*1./float(nt)
        print ' * Calendar: ', vtime[:]
    
        Xd_m_x_A = nmp.zeros((nt,nj,ni))
        for jt in range(nt):
            Xd_m_x_A[jt,:,:] = Xd_m[jt,:,:]*Xa[:,:]
            
        joce = 0
        for cocean in list_basin_names[:]:
            colnm = list_basin_lgnms[joce]
            print '   ===> '+cvar+' in basin '+cocean+' ('+colnm+')'
    
            # Decrasing the domain size if possible:
            (vjj , vji)  = nmp.where(mask[joce,:,:]>0.5)
            j1 = max( nmp.min(vjj)-2 , 0    )
            i1 = max( nmp.min(vji)-2 , 0    )
            j2 = min( nmp.max(vjj)+2 , nj-1 ) + 1
            i2 = min( nmp.max(vji)+2 , ni-1 ) + 1
            
            if (i1,i2) == (0,ni): i2 = i2-2 ; # Mind east-west periodicity overlap of 2 points...
            if (i1,j1,i2,j2) != (0,0,ni,nj): print '       ===> zooming on i1,j1 -> i2,j2:', i1,j1,'->',i2,j2
    
            v_x_m = nmp.zeros(nt)
            for jt in range(nt):
                v_x_m[jt] = nmp.sum( rmult*Xd_m_x_A[jt,j1:j2,i1:i2]*mask[joce,j1:j2,i1:i2] )
    
            cf_out = vdic['DIAG_D']+'/Sflx_'+cvar+'_'+CONFEXP+'_'+cocean+'.nc'
            bnc.wrt_appnd_1d_series(vtime, v_x_m, cf_out, cvar,
                                    cu_t='year', cu_d=cunit,
                                    cln_d ='2D-integral of '+cvar+' on region '+colnm)
            joce = joce + 1

        print '\n +++ '+cnexec+' => Done with 2D-summing of variable '+cvar+'!\n'
        del Xd_m



print '\n *** EXITING '+cnexec+' for year '+cyear+' and variable '+cvar+'!\n'
