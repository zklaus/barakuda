#!/usr/bin/env python

#       B a r a K u d a
#
#  Prepare 2D maps (monthly) that will later become a GIF animation!
#  NEMO output and observations needed
#
#    L. Brodeau, november 2016

import sys
import os
from string import replace
import numpy as nmp

from netCDF4 import Dataset

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import datetime

import barakuda_colmap as bcm

import barakuda_tool as bt

year_ref_ini = 1990

#jt0 = 248
jt0 = 0

#CTATM = 'T255'
CTATM = 'T1279'

cbox = 'NAtl'

fig_type='png'

narg = len(sys.argv)
if narg < 4: print 'Usage: '+sys.argv[0]+' <file> <variable> <LSM_file>'; sys.exit(0)
cf_in = sys.argv[1] ; cv_in=sys.argv[2] ; cf_lsm=sys.argv[3]

# Ice:
cv_ice  = 'siconc'
cf_ice = replace(cf_in, 'grid_T', 'icemod')
rmin_ice = 0.5
cpal_ice = 'ncview_bw'
vcont_ice = nmp.arange(rmin_ice, 1.05, 0.05)

if cv_in == 'sosstsst':
    cfield = 'SST'
    tmin=-2. ;  tmax=28.   ;  dtemp = 1.
    cpal_fld = 'ncview_nrl'    
    cunit = r'$^{\circ}C$'
    cb_jump = 2
    
elif cv_in == 'somxl010':
    cfield == 'MLD'
    tmin=50. ;  tmax=1500. ;  dtemp = 50.
    cpal_fld = 'viridis_r'
    

bt.chck4f(cf_ice)

bt.chck4f(cf_in)
id_fld = Dataset(cf_in)
vtime = id_fld.variables['time_counter'][:]
id_fld.close()
Nt = len(vtime)

bt.chck4f(cf_lsm)
id_lsm = Dataset(cf_lsm)
XMSK  = id_lsm.variables['tmask'][0,0,:,:] ; # t, y, x
id_lsm.close()
[ nj , ni ] = nmp.shape(XMSK)

pmsk = nmp.ma.masked_where(XMSK[:,:] > 0.2, XMSK[:,:]*0.+40.)





idx_oce = nmp.where(XMSK[:,:] > 0.5)






params = { 'font.family':'Ubuntu',
           'font.size':       int(12),
           'legend.fontsize': int(12),
           'xtick.labelsize': int(12),
           'ytick.labelsize': int(12),
           'axes.labelsize':  int(12) }
mpl.rcParams.update(params)
cfont_clb   = { 'fontname':'Arial', 'fontweight':'normal', 'fontsize':13 }
cfont_title = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':18 }
cfont_mail  = { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':9, 'color':'0.5' }


# Colormaps for fields:
pal_fld = bcm.chose_colmap(cpal_fld)
norm_fld = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)

pal_ice = bcm.chose_colmap(cpal_ice)
norm_ice = colors.Normalize(vmin = rmin_ice, vmax = 1, clip = False)

pal_lsm = bcm.chose_colmap('blk')
norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)

rh = 7.5

for jt in range(jt0,Nt):

    ct = '%3.3i'%(jt+1)

    cd = str(datetime.datetime.strptime('1990 '+ct, '%Y %j'))
    cdate = cd[:10] ; print ' *** cdate :', cdate

    cfig = 'figs/zoom_'+cv_in+'_NEMO'+'_d'+ct+'.'+fig_type    

    fig = plt.figure(num = 1, figsize=(rh,rh), dpi=None, facecolor='w', edgecolor='k')
    ax  = plt.axes([0.065, 0.05, 0.9, 1.], axisbg = 'k')

    vc_fld = nmp.arange(tmin, tmax + dtemp, dtemp)


    print "Reading record #"+str(ct)+" of "+cv_in+" in "+cf_in
    id_fld = Dataset(cf_in)
    XFLD  = id_fld.variables[cv_in][jt,:,:] ; # t, y, x
    id_fld.close()
    print "Done!"

    print "Ploting"
    cf = plt.imshow(XFLD[:,:], cmap = pal_fld, norm = norm_fld)
    del XFLD
    print "Done!"
    
    # Ice
    if not cfield == 'MLD':
        print "Reading record #"+str(ct)+" of "+cv_ice+" in "+cf_ice
        id_ice = Dataset(cf_ice)
        XICE  = id_ice.variables[cv_ice][jt,:,:] ; # t, y, x
        id_ice.close()
        print "Done!"

        #XM[:,:] = XMSK[:,:] 
        #bt.drown(XICE, XM, k_ew=2, nb_max_inc=10, nb_smooth=10)
        #ci = plt.contourf(XICE[:,:], vcont_ice, cmap = pal_ice, norm = norm_ice) #

        pice = nmp.ma.masked_where(XICE < rmin_ice, XICE)
        ci = plt.imshow(pice, cmap = pal_ice, norm = norm_ice) ; del pice, ci
        del XICE


    cm = plt.imshow(pmsk, cmap = pal_lsm, norm = norm_lsm)
    
    plt.axis([ 0, ni, 0, nj])

    plt.title('NEMO: '+cfield+', coupled ORCA12-'+CTATM+', '+cdate, **cfont_title)



    ax2 = plt.axes([0.055, 0.067, 0.93, 0.025])
    clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend='both')
    #clb = plt.colorbar(cf, ticks=vc_fld, orientation='horizontal', drawedges=False, pad=0.07, shrink=1., aspect=40) #
    cb_labs = [] ; cpt = 0
    for rr in vc_fld:
        if cpt % cb_jump == 0:
            cb_labs.append(str(int(rr)))
        else:
            cb_labs.append(' ')
        cpt = cpt + 1
    clb.ax.set_xticklabels(cb_labs)
    clb.set_label(cunit, **cfont_clb)

    del cf


    ax.annotate('laurent.brodeau@bsc.es', xy=(1, 4), xytext=(890, -150), **cfont_mail)
    
    plt.savefig(cfig, dpi=160, orientation='portrait', transparent=False)
    print cfig+' created!\n'
    plt.close(1)


    del cm, fig, ax, clb

