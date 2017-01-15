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

lsst=True
lshf=True


#venv_needed = {'ORCA','RUN','DIAG_D','MM_FILE','NN_SST','NN_T','NN_S','NN_ICEF',
#               'F_T_CLIM_3D_12','F_S_CLIM_3D_12','SST_CLIM_12','NN_SST_CLIM','NN_T_CLIM','NN_S_CLIM'}

#vdic = bt.check_env_var(sys.argv[0], venv_needed)

#CONFRUN = vdic['ORCA']+'-'+vdic['RUN']


fig_type='png'


if lsst:
    tmin=-20.  ;  tmax=12. ;  dt = 1.
    cf_sst  = '/home/laurent/tmp/T2M_ICMGG_C120_1990.nc4'
    csst  = 'T2M'
    cpal_sst = 'sstnw'
    cfield = 'SST'
    
if lshf:
    tmin=-1200. ;  tmax=200. ;  dt = 25.
    cf_sst   = '/home/laurent/tmp/IFS/C120_1990_1d_SNHF.nc4'
    csst     = 'SNHF'
    #cpal_sst = 'spectral'
    #cpal_sst = 'rainbow'
    cpal_sst = 'gist_ncar'
    #cpal_sst = 'nipy_spectral'
    cfield = 'Net Heat Flux'




cf_msk = '/home/laurent/tmp/LSM_ICMGG_T255.nc4'
cmsk  = 'LSM'

# South Greenland:
i1 = 412; i2 =486
j1 = 22 ; j2 = 56

#Global T255:
#i1 = 0 ; i2 =512
#j1 = 0 ; j2 = 256


bt.chck4f(cf_msk)
id_msk = Dataset(cf_msk)
XMSK  = id_msk.variables[cmsk][0,j1:j2,i1:i2] ; # t, y, x
id_msk.close()


[ nj , ni ] = nmp.shape(XMSK)


bt.chck4f(cf_sst)
id_sst = Dataset(cf_sst)
XSST   = id_sst.variables[csst][:,j1:j2,i1:i2]
id_sst.close()
[ Nt, nj0 , ni0 ] = nmp.shape(XSST)


if lsst: XSST[:,:,:] = XSST[:,:,:] - 273.15

#idx_oce = nmp.where(XMSK[:,:] > 0.01)

vc_sst = nmp.arange(tmin, tmax + dt, dt)

psst = nmp.zeros((nj,ni))


params = { 'font.family':'Ubuntu',
           'font.size':       int(15),
           'legend.fontsize': int(15),
           'xtick.labelsize': int(15),
           'ytick.labelsize': int(15),
           'axes.labelsize':  int(15) }
mpl.rcParams.update(params)

for jt in range(Nt):


    ct = '%3.3i'%(jt+1)

    cd = str(datetime.datetime.strptime('1990 '+ct, '%Y %j'))
    cdate = cd[:10] ; print ' *** cdate :', cdate

    cfig = csst+'_IFS'+'_d'+ct+'.'+fig_type
    
    #psst = nmp.ma.masked_where(XMSK[:,:] < 0.2, XSST[jt,:,:])
    #pice = nmp.ma.masked_where(XMSK[:,:] < 0.2, XICE[jt,:,:])

    psst[:,:] = XSST[jt,:,:]

    fig = plt.figure(num = 1, figsize=(10,8.5), dpi=None, facecolor='w', edgecolor='k')
    ax  = plt.axes([0.05, -0.05, 0.9, 0.999], axisbg = 'k')


    # Pal_Sst:
    pal_sst = bcm.chose_palette(cpal_sst)
    norm_sst = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)

    pal_msk = bcm.chose_palette('blk')
    norm_msk = colors.Normalize(vmin = 0., vmax = 1., clip = False)

    cf = plt.pcolor(nmp.flipud(psst), cmap = pal_sst, norm = norm_sst)


    # Mask
    pmsk = nmp.ma.masked_where(XMSK[:,:] < 0.5, XMSK[:,:]*0.+40.)
    plt.pcolor(nmp.flipud(pmsk), cmap = pal_msk, norm = norm_msk)
    
    plt.axis([ 0, ni, 0, nj])

    clb = plt.colorbar(cf, ticks=vc_sst, orientation='horizontal', drawedges=False, pad=0.07, shrink=1., aspect=40)
    cb_labs = [] ; cpt = 0
    for rr in vc_sst:
        if cpt % 4 == 0:
            cb_labs.append(str(int(rr)))
        else:
            cb_labs.append(' ')
        cpt = cpt + 1
    clb.ax.set_xticklabels(cb_labs)
    cfont = { 'fontname':'Arial', 'fontweight':'normal', 'fontsize':16 }
    clb.set_label(r'$W/m^2$', **cfont)

    cfont = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':24 }
    plt.title('IFS: '+cfield+', coupled ORCA12-T255, '+cdate)

    
    plt.savefig(cfig, dpi=120, orientation='portrait', transparent=False)
    print cfig+' created!\n'
    plt.close(1)
                




