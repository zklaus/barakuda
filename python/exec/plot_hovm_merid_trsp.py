#!/usr/bin/env python
#
# L. Brodeau, July 2011
#

import sys
import os
import numpy as nmp
from netCDF4 import Dataset

import barakuda_orca as bo
import barakuda_plot as bp
import barakuda_tool as bt


dy = 10. ; # latitude increment in figure...

venv_needed = {'ORCA','EXP','DIAG_D','BM_FILE'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

CONFEXP = vdic['ORCA']+'-'+vdic['EXP']


path_fig=vdic['DIAG_D']+'/'

fig_type='png'

cf_in = vdic['DIAG_D']+'/merid_transport_T_S_'+CONFEXP+'.nc'


if not os.path.exists(cf_in):
    print ' ERROR: plot_hovm_merid_trsp.py => old ascii file system not supported anymore!'
    print '        => need file '+cf_in+' !!!' ; sys.exit(0)



#list_basin_names, list_basin_lgnms = bo.get_basin_info(vdic['BM_FILE'])
# As in cdfmhst.F90:
list_basin_names = [ 'GLO','atl','pac','ind' ]
list_basin_lgnms = [ 'Global Ocean','Atlantic Ocean','Pacific Ocean','Indian Ocean' ]
nbasins = len(list_basin_names)


id_in = Dataset(cf_in)

vyear = id_in.variables['time'][:]   ; Nby = len(vyear) ; ittic = bt.iaxe_tick(Nby)
vyear = vyear - 0.5
vlat  = id_in.variables['lat'][:]    ; Nlat = len(vlat)

for joce in range(nbasins):

    cbas   = list_basin_names[joce] ; # name as in cf_in ...
    cbasin = list_basin_lgnms[joce] ; # long name of basin
    
    if joce == 0:
        Xheat = nmp.zeros(nbasins*Nby*Nlat) ; Xheat.shape = [ nbasins, Nby, Nlat ]
        Xsalt = nmp.zeros(nbasins*Nby*Nlat) ; Xsalt.shape = [ nbasins, Nby, Nlat ]
        rmiss_val = id_in.variables['zomht_'+cbas].getncattr('_FillValue')

    Xheat[joce,:,:] = id_in.variables['zomht_'+cbas][:,:]
    Xsalt[joce,:,:] = id_in.variables['zomst_'+cbas][:,:]
    print ' *** zomht_'+cbas+' and zomst_'+cbas+' sucessfully read into '+cf_in

id_in.close()

print ''




mask  = nmp.zeros((Nby, Nlat)) + 1.

for joce in range(nbasins):

    if joce == 0:
        vyear = nmp.trunc(vyear) + 0.5 ; # in case 1990 and not 1990.5 !!!
        yr1=float(int(min(vyear)))
        yr2=float(int(max(vyear)))

    cbas   = list_basin_names[joce] ; # name as in cf_in ...
    cbasin = list_basin_lgnms[joce] ; # long name of basin

    #mask[:,:] = 0
    #Lfinite = nmp.isfinite(Xheat[joce,:,:]) ; idx_good = nmp.where(Lfinite)
    #mask[idx_good] = 1
            
    # time record to consider to chose a min and a max for colorbar:
    jt_ini = 5
    if Nby <= 10: jt_ini = 1
    if Nby ==  1: jt_ini = 0

    # Finding the first and last j that make sense (not NaN)
    js = -1 ; lfound=False
    while not lfound:
        js = js + 1
        #if not nmp.isnan(Xheat[joce,1,js]): lfound=True
        if Xheat[joce,1,js] != rmiss_val:  lfound=True
    je = Nlat ; lfound=False
    while not lfound:
        je = je - 1
        #if not nmp.isnan(Xheat[joce,1,je]): lfound=True
        if Xheat[joce,1,je] != rmiss_val:  lfound=True
        
    dyh=dy/2
    ymin = (int(vlat[js]/dyh)-1)*dyh
    ymax = (int(vlat[je]/dyh)+1)*dyh

    # min and max for field:
    [ rmin, rmax, rdf ] = bt.get_min_max_df(Xheat[joce,jt_ini:,js+1:je-1],40)
        
    bp.plot("hovmoeller")(vyear[:], vlat[:], nmp.flipud(nmp.rot90(Xheat[joce,:,:])), mask,
                          rmin, rmax, rdf, c_y_is='latitude',
                          cpal='RdBu_r', tmin=yr1, tmax=yr2+1., dt=ittic, lkcont=True,
                          ymin = ymin, ymax = ymax, dy=dy,
                          cfignm=path_fig+'MHT_'+CONFEXP+'_'+cbas, cbunit='PW', ctunit='',
                          cyunit=r'Latitude ($^{\circ}$N)',
                          ctitle=CONFEXP+': Northward advective meridional heat transport, '+cbasin,
                          cfig_type=fig_type, i_cb_subsamp=2, l_y_increase=True)



    # Salt transport

    [ rmin, rmax, rdf ] = bt.get_min_max_df(Xsalt[joce,jt_ini:,js+1:je-1],40)

    bp.plot("hovmoeller")(vyear[:], vlat[:], nmp.flipud(nmp.rot90(Xsalt[joce,:,:])), mask,
                          rmin, rmax, rdf, c_y_is='latitude',
                          cpal='PiYG_r', tmin=yr1, tmax=yr2+1., dt=ittic, lkcont=True,
                          ymin = ymin, ymax = ymax, dy=dy,
                          cfignm=path_fig+'MST_'+CONFEXP+'_'+cbas, cbunit=r'10$^3$ tons/s', ctunit='',
                          cyunit=r'Latitude ($^{\circ}$N)',
                          ctitle=CONFEXP+': Northward advective meridional salt transport, '+cbasin,
                          cfig_type=fig_type, i_cb_subsamp=2, l_y_increase=True)
    

