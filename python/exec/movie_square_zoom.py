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


#venv_needed = {'ORCA','RUN','DIAG_D','MM_FILE','NN_SST','NN_T','NN_S','NN_ICEF',
#               'F_T_CLIM_3D_12','F_S_CLIM_3D_12','SST_CLIM_12','NN_SST_CLIM','NN_T_CLIM','NN_S_CLIM'}

#vdic = bt.check_env_var(sys.argv[0], venv_needed)

#CONFRUN = vdic['ORCA']+'-'+vdic['RUN']

tmin=-2. ;  tmax=12.   ;  dtemp = 1.
imin=0.  ;  imax=0.99  ;  dice = 0.1
cfield = 'SST'	

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


cf_sst  = '/home/laurent/tmp/SGL_C120_1d_19900101_19900731_grid_T.nc4'
csst  = 'sosstsst'

cf_ice  = '/home/laurent/tmp/SGL_C120_1d_19900101_19900731_icemod.nc4'
cice  = 'siconc'

cf_msk = '/home/laurent/tmp/ZOOMs/SGL_mesh_mask.nc4'
cmsk  = 'tmask'


#path_fig = 'movies'
 
#os.system("mkdir -p "+path_fig)


bt.chck4f(cf_msk)
id_msk = Dataset(cf_msk)
XMSK  = id_msk.variables[cmsk][0,0,:,:] ; # t, y, x
id_msk.close()


[ nj , ni ] = nmp.shape(XMSK)


cpal_sst = 'sstnw'
bt.chck4f(cf_sst)
id_sst = Dataset(cf_sst)
XSST  = id_sst.variables[csst][:,:,:] ; # t, y, x
id_sst.close()
[ Nt, nj0 , ni0 ] = nmp.shape(XSST)


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
           'ytick.labelsize': int(15),
           'axes.labelsize':  int(15) }
mpl.rcParams.update(params)

idx_oce = nmp.where(XMSK[:,:] > 0.5)


psst = nmp.zeros((nj,ni))
pice = nmp.zeros((nj,ni))

for jt in range(Nt):


    ct = '%3.3i'%(jt+1)

    cd = str(datetime.datetime.strptime('1990 '+ct, '%Y %j'))
    cdate = cd[:10] ; print ' *** cdate :', cdate

    cfig = csst+'_NEMO'+'_d'+ct+'.'+fig_type
    
    #psst = nmp.ma.masked_where(XMSK[:,:] < 0.2, XSST[jt,:,:])
    #pice = nmp.ma.masked_where(XMSK[:,:] < 0.2, XICE[jt,:,:])

    psst[:,:] = XSST[jt,:,:]

    pice[:,:] = XICE[jt,:,:]
    bt.drown(pice, XMSK, k_ew=2, nb_max_inc=10, nb_smooth=10)
    
    vc_sst = nmp.arange(tmin, tmax + dtemp, dtemp)

    fig = plt.figure(num = 1, figsize=(10,9), dpi=None, facecolor='w', edgecolor='k')
    ax  = plt.axes([0.05, -0.05, 0.9, 0.999], axisbg = 'k')


    # Pal_Sst:
    pal_sst = bcm.chose_palette(cpal_sst)
    norm_sst = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)
    pal_ice = bcm.chose_palette(cpal_ice)
    norm_ice = colors.Normalize(vmin = imin, vmax = imax, clip = False)

    pal_msk = bcm.chose_palette('blk')
    norm_msk = colors.Normalize(vmin = 0., vmax = 1., clip = False)

    

    cf = plt.pcolor(psst, cmap = pal_sst, norm = norm_sst)

    #plt.pcolor(pice, cmap = pal_ice, norm = norm_ice)
    plt.contourf(pice, [0.25,0.5,0.75,1.], cmap = pal_ice, norm = norm_ice)


    # Mask
    pmsk = nmp.ma.masked_where(XMSK[:,:] > 0.2, XMSK[:,:]*0.+40.)
    plt.pcolor(pmsk, cmap = pal_msk, norm = norm_msk)


    
    plt.axis([ 0, ni, 0, nj])

    clb = plt.colorbar(cf, ticks=vc_sst, orientation='horizontal', drawedges=False, pad=0.07, shrink=1., aspect=40)
    #cb_labs = [] ; cpt = 0
    #for rr in vc_sst:
    #    if cpt % 4 == 0:
    #        cb_labs.append(str(int(rr)))
    #    else:
    #        cb_labs.append(' ')
    #    cpt = cpt + 1
    #clb.ax.set_xticklabels(cb_labs)
    cfont = { 'fontname':'Arial', 'fontweight':'normal', 'fontsize':16 }
    clb.set_label(r'$^{\circ}C$', **cfont)

    cfont = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':24 }
    plt.title('NEMO: '+cfield+', coupled ORCA12-T255, '+cdate)

    
    plt.savefig(cfig, dpi=120, orientation='portrait', transparent=False)
    print cfig+' created!\n'
    plt.close(1)
                




