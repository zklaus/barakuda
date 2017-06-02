#!/usr/bin/env python

#       B a r a K u d a
#
#     Generate global plot of the barotropic stream function
#
#       L. Brodeau, 2017

import sys
import numpy as nmp
from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_plot as bp

cv_psi = 'sobarstf'

venv_needed = {'ORCA','EXP','DIAG_D','MM_FILE'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFEXP = vdic['ORCA']+'-'+vdic['EXP']

path_fig='./'
fig_type='png'



narg = len(sys.argv)
if narg < 3: print 'Usage: '+sys.argv[0]+' <year1> <year2>'; sys.exit(0)
cy1 = sys.argv[1] ; cy2=sys.argv[2]; jy1=int(cy1); jy2=int(cy2)


jy1_clim = jy1 ; jy2_clim = jy2
print ' => mean on the clim : ', jy1_clim, jy2_clim, '\n'


# Getting coordinates:
bt.chck4f(vdic['MM_FILE'])
id_mm = Dataset(vdic['MM_FILE'])
xlon   = id_mm.variables['glamt'][0,:,:] ; xlat = id_mm.variables['gphit'][0,:,:]
Xmask = id_mm.variables['tmask'][0,0,:,:]
#Xe1t = id_mm.variables['e1t'][0,:,:]
#Xe2t = id_mm.variables['e2t'][0,:,:]
id_mm.close()


#  Getting NEMO mean monthly climatology of PSI:
cf_nemo_mnmc = vdic['DIAG_D']+'/clim/mclim_'+CONFEXP+'_'+cy1+'-'+cy2+'_PSI.nc4'

bt.chck4f(cf_nemo_mnmc)
id_nemo = Dataset(cf_nemo_mnmc)
psi   = 1.E-6 * id_nemo.variables[cv_psi][:,:,:]
id_nemo.close()

[ nt, nj, ni ] = psi.shape ; print ' Shape of Psi :', nt, nj, ni, '\n'

psi_plot = nmp.zeros((nj,ni))


if nt == 1:
    psi_plot[:,:] = psi[0,:,:]
elif nt > 1:
    psi_plot[:,:] = nmp.mean(psi[:,:,:],axis=0)
else:
    print ' psi.py : Problem!!!'


print psi_plot[61,:]


#ztot = nmp.sum(psi_plot*Xmask*Xe1t*Xe2t)/nmp.sum(Xmask*Xe1t*Xe2t)
#print 'ztot =', ztot
#
#psi_plot = psi_plot - ztot
#cztot = str(round(ztot,2))


# the Jean-Marc Molines method:
ji_lat0 = nmp.argmax(xlat[nj-1,:])

bp.plot("2d")(xlon[0,:], xlat[:,ji_lat0], psi_plot[:,:], Xmask, -100., 100., 5.,
              corca=vdic['ORCA'], lkcont=True, cpal='BrBG_r',
              cfignm=path_fig+'psi_mean_'+CONFEXP, cbunit=r'$(10^{6} m^3/s)$',
              ctitle='Mean barotropic stream function , '+CONFEXP+' ('+cy1+'-'+cy2+')',
              lforce_lim=True, i_cb_subsamp=2,
              cfig_type=fig_type, lat_min=-70., lat_max=68., lpix=False, vcont_spec = [ 0. ])

