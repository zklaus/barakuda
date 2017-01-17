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

from netCDF4 import Dataset

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import datetime

import barakuda_colmap as bcm

import barakuda_tool as bt
import barakuda_plot as bp


imin=0.  ;  imax=0.99  ;  dice = 0.1

cfield = 'SST' ; cfld  = 'sosstsst' ; cpal_fld = 'sstnw'
cfield = 'MLD' ; cfld  = 'somxl010' ; cpal_fld = 'viridis_r'

fig_type='png'



#narg = len(sys.argv)
#if narg < 4:
#    print 'Usage: '+sys.argv[0]+' <NEMO file (1 year, monthyly)> <year>'
#    print '          ...'
#    sys.exit(0)

#cf_in = sys.argv[1]
#cy    = sys.argv[2] ; jy=int(cy)
#cvar  = sys.argv[3]

#if not cvar in ['sst','sss','ice']:
#    print 'ERROR (prepare_movies.py): variable '+cvar+' not supported yet!'
#    sys.exit(0)

#cbox  = 'SGL'
#croot = cbox+'_C120_1d_19900101_19900930'
cbox = 'NAtl'
croot = cbox+'_C120_1d_19901101_19901130'

cf_fld  = '/home/laurent/tmp/NEMO/'+croot+'_grid_T.nc4'


cf_ice  = '/home/laurent/tmp/NEMO/'+croot+'_icemod.nc4'
cice  = 'siconc'

cf_msk = '/home/laurent/tmp/NEMO/ZOOMs/'+cbox+'_mesh_mask.nc4'
cmsk  = 'tmask'


#path_fig = 'movies'
 
#os.system("mkdir -p "+path_fig)


bt.chck4f(cf_msk)
id_msk = Dataset(cf_msk)
XMSK  = id_msk.variables[cmsk][0,0,:,:] ; # t, y, x
id_msk.close()


[ nj , ni ] = nmp.shape(XMSK)



bt.chck4f(cf_fld)
id_fld = Dataset(cf_fld)
XFLD  = id_fld.variables[cfld][:,:,:] ; # t, y, x
id_fld.close()
[ Nt, nj0 , ni0 ] = nmp.shape(XFLD)


if not cfield == 'MLD':
    cpal_ice = 'ice'
    bt.chck4f(cf_ice)
    id_ice = Dataset(cf_ice)
    XICE  = id_ice.variables[cice][:,:,:] ; # t, y, x
    id_ice.close()
    [ Nt, nj0 , ni0 ] = nmp.shape(XICE)



params = { 'font.family':'Ubuntu',
           'font.size':       int(15),
           'legend.fontsize': int(15),
           'xtick.labelsize': int(15),
           'ytick.labelsize': int(15), # 
           'axes.labelsize':  int(15) }
mpl.rcParams.update(params)

idx_oce = nmp.where(XMSK[:,:] > 0.5)


if not cfield == 'MLD': pice = nmp.zeros((nj,ni))

for jt in range(Nt):


    ct = '%3.3i'%(jt+1)

    cd = str(datetime.datetime.strptime('1990 '+ct, '%Y %j'))
    cdate = cd[:10] ; print ' *** cdate :', cdate

    cfig = 'figs/'+cfld+'_NEMO'+'_d'+ct+'.'+fig_type    

    if cbox == 'SGL':
        if cfield == 'SST': tmin=-2. ;  tmax=12.   ;  dtemp = 1.
        fig = plt.figure(num = 1, figsize=(10,9), dpi=None, facecolor='w', edgecolor='k')
        ax  = plt.axes([0.05, -0.06, 0.93, 1.02], axisbg = 'k')
    elif cbox == 'NAtl':
        if cfield == 'SST': tmin=-2. ;  tmax=26.   ;  dtemp = 1.
        if cfield == 'MLD': tmin=50. ;  tmax=1500. ;  dtemp = 50.
        fig = plt.figure(num = 1, figsize=(10,10), dpi=None, facecolor='w', edgecolor='k')
        ax  = plt.axes([0.051, -0.06, 0.92, 1.02], axisbg = 'k')

    vc_fld = nmp.arange(tmin, tmax + dtemp, dtemp)

    # Pal_fld:
    pal_fld = bcm.chose_palette(cpal_fld)
    norm_fld = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)

    
    cf = plt.pcolor(XFLD[jt,:,:], cmap = pal_fld, norm = norm_fld)

    
    # Ice
    if not cfield == 'MLD':
        pal_ice = bcm.chose_palette(cpal_ice)
        norm_ice = colors.Normalize(vmin = imin, vmax = imax, clip = False)
        pice[:,:] = XICE[jt,:,:]
        bt.drown(pice, XMSK, k_ew=2, nb_max_inc=10, nb_smooth=10)
        plt.contourf(pice, [0.25,0.5,0.75,1.], cmap = pal_ice, norm = norm_ice) # 


    # Mask
    pal_msk = bcm.chose_palette('blk')
    norm_msk = colors.Normalize(vmin = 0., vmax = 1., clip = False)
    
    pmsk = nmp.ma.masked_where(XMSK[:,:] > 0.2, XMSK[:,:]*0.+40.)
    plt.pcolor(pmsk, cmap = pal_msk, norm = norm_msk)


    
    plt.axis([ 0, ni, 0, nj])

    clb = plt.colorbar(cf, ticks=vc_fld, orientation='horizontal', drawedges=False, pad=0.07, shrink=1., aspect=40) # 

    if cfield == 'MLD':         # 
        cb_labs = [] ; cpt = 0
        for rr in vc_fld:
            if (cpt+1) % 2 == 0:
                cb_labs.append(str(int(rr)))
            else:
                cb_labs.append(' ')
            cpt = cpt + 1
        clb.ax.set_xticklabels(cb_labs)
        
    cfont = { 'fontname':'Arial', 'fontweight':'normal', 'fontsize':16 }
    clb.set_label(r'$^{\circ}C$', **cfont)

    cfont = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':22 }
    plt.title('NEMO: '+cfield+', coupled ORCA12-T255, '+cdate, **cfont)

    
    plt.savefig(cfig, dpi=120, orientation='portrait', transparent=False)
    print cfig+' created!\n'
    plt.close(1)
                




