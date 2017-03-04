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


venv_needed = {'ORCA','EXP','DIAG_D','DENSITY_SECTION_FILE'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFEXP = vdic['ORCA']+'-'+vdic['EXP']

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

vtime  = id_in.variables['time'][:]
(nby, nbm, nbr, ittic) = bt.test_nb_years(vtime, 'trsp_sigma')

vsigma = id_in.variables['sigma_bins'][:] ; nbin = len(vsigma)
print '      => '+str(nbin)+' sigma-density bins and '+str(nbr)+' time records...'



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

    # Annual array: lolo

    if nbm >= 12:
        # the file contains monthly data (nbm=-1 otherwize)
        vtime_ann, Xst_ann = bt.monthly_2_annual(vtime, Xst)
    else:
        vtime_ann = vtime
        Xst_ann   = Xst



    # FIGURE 1
    ###########
    rmax = nmp.amax(nmp.abs(Xst_ann)) ; rmin = -rmax
    rmin, rmax, dc = bp.__suitable_axis_dx__(rmin, rmax, nb_val=25, lsym0=True)
    i_cbssmp=2
    if dc in [0.5,1.,2.,5.]: i_cbssmp=1

    bp.plot("trsp_sig_class")(vtime_ann, vsigma_bounds, Xst_ann, rmin, rmax, dc,
                              lkcont=False, cpal='ncview_jaisnb', dt=ittic,
                              cfignm='transport_sigma_class_'+csec+'_'+CONFEXP,
                              cfig_type='png', ctitle=r'Transport by $\sigma_0$ class, '+csec+', '+CONFEXP,
                              vcont_spec1 = [], i_cb_subsamp=i_cbssmp)
    

    # Volume transport for density > rsigma_dense
    # ====================================

    if jsec == 0:
        v278 = nmp.zeros(nby*nbsec) ; v278.shape = [ nbsec, nby ]
        j278 = 0 ; nsig = len(vsigma)
        while vsigma[j278] < rsigdense0: j278 = j278 + 1
        
    for jt in range(nby): v278[jsec,jt] = nmp.sum(Xst_ann[j278:,jt])


    jsec = jsec + 1


# Closing netcdf file:
id_in.close()


bp.plot("1d_multi")(vtime_ann, v278, list_sections, cfignm='tr_sigma_gt278_'+CONFEXP,
                    dt=ittic, cyunit='Sv',
                    ctitle=r'Transport of volume for $\sigma_0$ > '+str(rsigdense0)+', '+CONFEXP,
                    ymin=0., ymax=0., loc_legend='out')


print '\n trsp_sigma.py => done!\n\n'

