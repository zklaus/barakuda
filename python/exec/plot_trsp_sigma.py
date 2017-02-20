#!/usr/bin/env python

#       B a r a K u d a
#
#     Generate time-series of volume trasnport by sigma0 class
#
#       L. Brodeau, 2013
#

import sys
import numpy as nmp

from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_plot as bp
import barakuda_physics as bphy


venv_needed = {'ORCA','RUN','DIAG_D','DENSITY_SECTION_FILE'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFRUN = vdic['ORCA']+'-'+vdic['RUN']

rsigdense0 = bphy.rsigma_dense

cfig_type = 'png'

path_fig =  vdic['DIAG_D']+'/'

cf_dens_sect =  vdic['DENSITY_SECTION_FILE']
print '  Using cf_dens_sect = '+cf_dens_sect
list_sections = bt.get_sections_from_file(cf_dens_sect)
print 'List of sections to treat: ', list_sections
nbsec = len(list_sections)



cf_in =  vdic['DIAG_D']+'/transport_by_sigma_class.nc' ; bt.chck4f(cf_in)



id_in = Dataset(cf_in)

vtime  = id_in.variables['time'][:]       ; nbm  = len(vtime)
vsigma = id_in.variables['sigma_bins'][:] ; nbin = len(vsigma)

print '      => '+str(nbin)+' sigma-density bins and '+str(nbm)+' time snapshots...'

if nbm % 12 != 0: print 'nbm is not a multiple of 12!'; sys.exit(0)


# Reconstructing bounds of bins:
vsigma_bounds = nmp.zeros(nbin+1)
dsigma = vsigma[1]-vsigma[0]
vsigma_bounds[:nbin] = vsigma[:] - 0.5*dsigma
vsigma_bounds[nbin]  = vsigma[nbin-1] + 0.5*dsigma



# Loop along sections:
######################

jsec = 0

for csec in list_sections:


    Xst = nmp.flipud(nmp.rot90(id_in.variables['sigtrsp_'+csec][:,:]))
    print ' Shape of "sigtrsp_'+csec+'" => ', nmp.shape(Xst)

    # Annual array:
    vtime_ann, Xst_ann = bt.monthly_2_annual(vtime, Xst)



    # FIGURE 1
    ###########

    ittic = bt.iaxe_tick(nbm/12)

    # We want rmax to be a multile of 0.2:
    rmax = nmp.amax(nmp.abs(Xst_ann))
    r1   = round(rmax+0.05,1)*100; rmax = (r1 + r1%20)/100. ; rmin = -rmax
    dc = max( (int(round(100*(rmax+0.6)/20.,2))/5*5)/100. , 0.01 )    
    
    bp.plot("trsp_sig_class")(vtime_ann, vsigma_bounds, Xst_ann, rmin, rmax, dc,
                              lkcont=False, cpal='bbr2_r', dt=ittic,
                              cfignm='transport_sigma_class_'+csec+'_'+CONFRUN,
                              cfig_type='png', ctitle=r'Transport by $\sigma_0$ class, '+csec+', '+CONFRUN,
                              vcont_spec1 = [], i_cb_subsamp=2)
    

    # Volume transport for density > rsigma_dense
    # ====================================

    if jsec == 0:
        v278 = nmp.zeros(nbm/12*nbsec) ; v278.shape = [ nbsec, nbm/12 ]
        j278 = 0 ; nsig = len(vsigma)
        while vsigma[j278] < rsigdense0: j278 = j278 + 1
        
    for jt in range(nbm/12): v278[jsec,jt] = nmp.sum(Xst_ann[j278:,jt])


    jsec = jsec + 1


# Closing netcdf file:
id_in.close()



ittic = bt.iaxe_tick(nbm/12)

bp.plot("1d_multi")(vtime_ann, v278, list_sections, cfignm='tr_sigma_gt278_'+CONFRUN,
                    dt=ittic, cyunit='Sv',
                    ctitle=r'Transport of volume for $\sigma_0$ > '+str(rsigdense0)+', '+CONFRUN,
                    ymin=0., ymax=0.)


print '\n trsp_sigma.py => done!\n\n'

