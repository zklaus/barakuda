#!/usr/bin/env python

#       B a r a K u d a
#
#  Prepare 2D maps (monthly) that will later become a GIF animation!
#  NEMO output and observations needed
#
#    L. Brodeau, May 2017

import sys
import os
from string import replace
import numpy as nmp

from netCDF4 import Dataset

#from mpl_toolkits.basemap import Basemap
#from mpl_toolkits.basemap import shiftgrid
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import warnings
warnings.filterwarnings("ignore")

import datetime

import barakuda_colmap as bcm

import barakuda_tool as bt


CNEMO = 'NATL60'

#color_top = 'white'
color_top = 'k'

l_do_ice = False

year_ref_ini = 2013

#jt0 = 248
jt0 = 0




i2=0
j2=0

#cbox = 'Biscay' ; i1 = 2880 ; j1 = 1060 ; rfact_zoom = 1. ; ny_res = 1080 ; nx_res = 1920
#cbox = 'NAtlt' ; i1 = 0 ; j1 = 150 ; rfact_zoom = 0.36  ; ny_res = 1080 ; nx_res = 1920

cbox = 'Med+NorthSea' ; i1 = 4086 ; j1 = 603 ; i2 = 5421 ; j2 = 3119 ; rfact_zoom = 0.3319

if i2>0 and j2>0:
    nx_res = i2-i1
    ny_res = j2-j1
else:
    i2 = i1 + nxr
    j2 = j1 + nyr

print ' i1,i2,j1,j2 =>', i1,i2,j1,j2

yx_ratio = float(ny_res)/float(nx_res)
nxr = int(1./rfact_zoom*nx_res)
nyr = int(1./rfact_zoom*ny_res)




fig_type='png'

narg = len(sys.argv)
if narg < 5: print 'Usage: '+sys.argv[0]+' <file> <variable> <LSM_file> <month>'; sys.exit(0)
cf_in = sys.argv[1] ; cv_in=sys.argv[2] ; cf_lsm=sys.argv[3] ; cmn=sys.argv[4]

# Ice:
if l_do_ice:
    cv_ice  = 'siconc'
    cf_ice = replace(cf_in, 'grid_T', 'icemod')
    rmin_ice = 0.5
    cpal_ice = 'ncview_bw'
    vcont_ice = nmp.arange(rmin_ice, 1.05, 0.05)

if cv_in == 'sosstsst':
    cfield = 'SST'
    tmin=0. ;  tmax=25.   ;  dtemp = 1.
    cpal_fld = 'ncview_nrl'    
    cunit = r'SST ($^{\circ}$C)'
    cb_jump = 1

if cv_in == 'sossheig':
    cfield = 'SSH'
    tmin=-0.3 ;  tmax=0.3   ;  dtemp = 0.1
    cpal_fld = 'ncview_jaisnc'    
    cunit = r'SSH (m)'
    cb_jump = 1

elif cv_in == 'somxl010':
    cfield == 'MLD'
    tmin=50. ;  tmax=1500. ;  dtemp = 50.
    cpal_fld = 'viridis_r'
    

if l_do_ice: bt.chck4f(cf_ice)

bt.chck4f(cf_in)
id_fld = Dataset(cf_in)
vtime = id_fld.variables['time_counter'][:]
id_fld.close()
Nt = len(vtime)

bt.chck4f(cf_lsm)
id_lsm = Dataset(cf_lsm)
XMSK  = id_lsm.variables['tmask'][0,0,j1:j2,i1:i2] ; # t, y, x
id_lsm.close()
[ nj , ni ] = nmp.shape(XMSK)

pmsk = nmp.ma.masked_where(XMSK[:,:] > 0.2, XMSK[:,:]*0.+40.)





idx_oce = nmp.where(XMSK[:,:] > 0.5)


#params = { 'font.family':'Ubuntu',
params = { 'font.family':'Helvetica Neue',
           'font.weight':    'normal',
           'font.size':       int(22*rfact_zoom**0.2),
           'legend.fontsize': int(22*rfact_zoom**0.2),
           'xtick.labelsize': int(22*rfact_zoom**0.2),
           'ytick.labelsize': int(22*rfact_zoom**0.2),
           'axes.labelsize':  int(22*rfact_zoom**0.2) }
mpl.rcParams.update(params)
cfont_clb  = { 'fontname':'Helvetica Neue', 'fontweight':'medium', 'fontsize':int(22*rfact_zoom**0.2), 'color':color_top }
cfont_date = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(22*rfact_zoom**0.2), 'color':'w' }
cfont_mail = { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14*rfact_zoom**0.2), 'color':'0.8'}
cfont_titl = { 'fontname':'Helvetica Neue', 'fontweight':'light', 'fontsize':int(80*rfact_zoom**0.2), 'color':'w' }


# Colormaps for fields:
pal_fld = bcm.chose_colmap(cpal_fld)
norm_fld = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)

if l_do_ice:
    pal_ice = bcm.chose_colmap(cpal_ice)
    norm_ice = colors.Normalize(vmin = rmin_ice, vmax = 1, clip = False)

pal_lsm = bcm.chose_colmap('land_dark')
norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)

rh = 16.*rfact_zoom

for jt in range(jt0,Nt):

    ch = '%2.2i'%((jt+1)%24)
    cd = '%3.3i'%(1+(jt+1)/24)
    
    ct = str(datetime.datetime.strptime(str(year_ref_ini)+'-'+cmn+'-'+cd+' '+ch, '%Y-%m-%j %H'))
    #print ' cmn = ', cmn
    #print ' ct = ', ct
    ct=ct[:5]+cmn+ct[7:] #lolo bug !!! need to do that to get the month and not "01"
    print ' ct = ', ct
    cday  = ct[:10]   ; print ' *** cday  :', cday
    chour = ct[11:13] ; print ' *** chour :', chour



    cfig = 'figs/zoom_'+cv_in+'_NEMO'+'_'+cday+'_'+chour+'.'+fig_type    

    #fig = plt.figure(num = 1, figsize=(rh,rh*(9./16.)), dpi=None, facecolor='w', edgecolor='0.5')
    fig = plt.figure(num = 1, figsize=(rh,rh*yx_ratio), dpi=None, facecolor='w', edgecolor='0.5')

    #ax  = plt.axes([0.065, 0.05, 0.9, 1.], axisbg = '0.5')
    ax  = plt.axes([0., 0., 1., 1.], axisbg = '0.5')

    vc_fld = nmp.arange(tmin, tmax + dtemp, dtemp)


    print "Reading record #"+str(jt)+" of "+cv_in+" in "+cf_in
    id_fld = Dataset(cf_in)
    XFLD  = id_fld.variables[cv_in][jt,j1:j2,i1:i2] ; # t, y, x
    id_fld.close()
    print "Done!"

    print '  *** dimension of array => ', nmp.shape(XFLD)

    print "Ploting"
    cf = plt.imshow(XFLD[:,:], cmap = pal_fld, norm = norm_fld)
    del XFLD
    print "Done!"
    
    # Ice
    if not cfield == 'MLD' and l_do_ice:
        print "Reading record #"+str(jt)+" of "+cv_ice+" in "+cf_ice
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

    #plt.title('NEMO: '+cfield+', coupled '+CNEMO+', '+cday+' '+chour+':00', **cfont_title)



    #ax2 = plt.axes([0.3, 0.08, 0.4, 0.025])

    if yx_ratio > 1.:
        ax2 = plt.axes([0.01, 0.08, 0.98, 0.025])
    else:
        ax2 = plt.axes([0.22, 0.08, 0.56, 0.025])
        
    clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend='both')
    if cb_jump > 1:
        cb_labs = [] ; cpt = 0
        for rr in vc_fld:
            if cpt % cb_jump == 0:
                if dtemp >= 1.: cb_labs.append(str(int(rr)))
                if dtemp <  1.: cb_labs.append(str(rr))
            else:
                cb_labs.append(' ')
            cpt = cpt + 1
        clb.ax.set_xticklabels(cb_labs)    
    clb.set_label(cunit, **cfont_clb)
    clb.ax.yaxis.set_tick_params(color=color_top) ; # set colorbar tick color    
    clb.outline.set_edgecolor(color_top) ; # set colorbar edgecolor         
    plt.setp(plt.getp(clb.ax.axes, 'xticklabels'), color=color_top) ; # set colorbar ticklabels
        
    del cf


    #x_annot = nxr-nxr*0.22*rfact_zoom**0.2 ; y_annot = 150
    x_annot = 650 ; y_annot = 1035

    ax.annotate('Date: '+cday+' '+chour+':00',   xy=(1, 4), xytext=(x_annot,    y_annot), **cfont_date)

    ax.annotate('laurent.brodeau@ocean-next.fr', xy=(1, 4), xytext=(x_annot+150, 20), **cfont_mail)

    ax.annotate('NATL60', xy=(1, 4), xytext=(int(0.025*nxr), int(0.7*nyr)), **cfont_titl)

    
    plt.savefig(cfig, dpi=120, orientation='portrait', facecolor='k')
    print cfig+' created!\n'
    plt.close(1)


    del cm, fig, ax, clb

