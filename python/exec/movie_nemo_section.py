#!/usr/bin/env python

#       B a r a K u d a
#
#  Prepare 2D maps (monthly) that will later become a GIF animation!
#  NEMO output and observations needed
#
#    L. Brodeau, May 2018

import sys
import os
from string import replace
import numpy as nmp

from netCDF4 import Dataset

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import warnings
warnings.filterwarnings("ignore")

from calendar import isleap
import datetime

import barakuda_colmap as bcm

import barakuda_tool as bt


vmn = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
vml = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]


color_top = 'white'
#color_top = 'k'



#jt0 = 248
jt0 = 0


j2=0
k2=0
l_show_cb = True
l_show_date = True
l_log_field = False
l_pow_field = False
l_annotate_name = False


narg = len(sys.argv)
if narg < 6: print 'Usage: '+sys.argv[0]+' <NEMOCONF> <file> <variable> <LSM_file> <YYYYMMDD (start)>'; sys.exit(0)
CNEMO  = sys.argv[1]
cf_in = sys.argv[2] ; cv_in=sys.argv[3]
cf_lsm=sys.argv[4] ; cf_date0=sys.argv[5]


if CNEMO == 'eNATL60':
    Nk0 = 300
    Nj0 = 4729-1
    l_show_cb = False
    l_show_date = True
    cdt = '1h'
    #cbox = 'FullMed' ; j1=5400 ; k1=1530 ; j2=Nj0 ; k2=3310 ; rfact_zoom = 0.79 ; vcb = [0.5, 0.875, 0.485, 0.02] ; font_rat = 2.*rfact_zoom ; l_annotate_name=False
    cbox = 'ALL' ; j1=0 ; k1=0 ; j2=Nj0 ; k2=Nk0 ; rfact_zoom = 0.3047 ; vcb = [0.59, 0.1, 0.38, 0.018] ; font_rat = 8.*rfact_zoom
    #cbox = 'Portrait' ; j1=2760 ; k1=1000 ; j2=4870 ; k2=4000 ; rfact_zoom = 1. ; vcb = [0.59, 0.1, 0.38, 0.018] ; font_rat = 1.*rfact_zoom ; l_annotate_name=False; l_show_date=False
    x_date = 1900 ; y_date = 20 ; # where to put the date

if CNEMO == 'NATL60':
    Nk0 = 300
    Nj0 = 3454-1
    #l_pow_field = True ; pow_field = 1.5
    l_show_cb = False
    l_show_date = False
    cdt = '1h'
    #cbox = 'zoom1' ; j1 = 1800 ; k1 = 950 ; j2 = j1+1920 ; k2 = k1+1080 ; rfact_zoom = 1. ; vcb = [0.5, 0.875, 0.485, 0.02] ; font_rat = 8.*rfact_zoom ; l_show_lsm = False
    #cbox = 'zoom1' ; j1 = 1800 ; k1 = 950 ; j2 = j1+2560 ; k2 = k1+1440 ; rfact_zoom = 1. ; vcb = [0.5, 0.875, 0.485, 0.02] ; font_rat = 8.*rfact_zoom
    cbox = 'ALL' ; j1=0 ; k1=0 ; j2=Nj0 ; k2=Nk0 ; rfact_zoom = 1. ; vcb = [0.59, 0.1, 0.38, 0.018] ; font_rat = 4.*rfact_zoom
    x_date = 350 ; y_date = 7 ; # where to put the date


if CNEMO == 'NANUK025':
    cdt = '3h'; cbox = 'ALL' ; j1 = 0 ; k1 = 0 ; j2 = 492 ; k2 = 614 ; rfact_zoom = 2. ; vcb = [0.5, 0.875, 0.485, 0.02] ; font_rat = 8.*rfact_zoom
    x_date = 350 ; y_date = 7 ; # where to put the date


print '\n rfact_zoom = ', rfact_zoom
print ' font_rat = ', font_rat, '\n'

nx_res = j2-j1
ny_res = k2-k1

print ' *** nx_res, ny_res =', nx_res, ny_res


print ' j1,j2,k1,k2 =>', j1,j2,k1,k2

yx_ratio = float(ny_res)/float(nx_res)

nxr = int(rfact_zoom*nx_res) ; # widt image (in pixels)
nyr = int(rfact_zoom*ny_res) ; # height image (in pixels)




dpi = 110

rh = round(float(nxr)/float(dpi),3) ; # width of figure as for figure...

fig_type='png'



cyr0=cf_date0[0:4]
cmn0=cf_date0[4:6]
cdd0=cf_date0[6:8]


l_3d_field = False


if cv_in in ['sosstsst','tos']:
    cfield = 'SST'
    tmin=0. ;  tmax=30.   ;  df = 2. ; cpal_fld = 'ncview_nrl'
    #tmin=0. ;  tmax=32.   ;  df = 2. ; cpal_fld = 'viridis'
    #tmin=4. ;  tmax=20.   ;  df = 1. ; cpal_fld = 'PuBu'
    cunit = r'SST ($^{\circ}$C)'
    cb_jump = 2
    l_show_cb = True

if cv_in == 'sossheig':
    cfield = 'SSH'
    #tmin=-0.5 ;  tmax=0.5   ;  df = 0.05
    #tmin=-1.2 ;  tmax=2.3   ;  df = 0.05 ; l_apply_lap = True
    #cpal_fld = 'ncview_jaisnc'
    #cpal_fld = 'PuBu'
    #cpal_fld = 'RdBu'
    #cpal_fld = 'BrBG'
    #
    #cpal_fld = 'on3' ; tmin=-1.2 ;  tmax=2.3   ;  df = 0.05 ; l_apply_lap = True
    cpal_fld = 'on2' ; tmin=-1.2 ;  tmax=1.2   ;  df = 0.05
    #cpal_fld = 'coolwarm' ; tmin=-1. ;  tmax=1.   ;  df = 0.05 ; l_apply_lap = True
    #cpal_fld = 'RdBu_r' ; tmin=-0.9 ;  tmax=-tmin   ;  df = 0.05 ; l_apply_lap = True
    #cpal_fld = 'gray_r' ; tmin=-0.3 ;  tmax=0.3   ;  df = 0.05 ; l_apply_lap = True
    #cpal_fld = 'bone_r' ; tmin=-0.9 ;  tmax=-tmin   ;  df = 0.05 ; l_apply_lap = True ; l_pow_field = True ; pow_field = 2.
    cunit = r'SSH (m)'
    cb_jump = 1

if cv_in == 'socurloverf':
    cfield = 'RV'
    cpal_fld = 'on2' ; tmin=-1. ;  tmax=1. ;  df = 0.05
    cunit = ''
    cb_jump = 1


if cv_in == 'vozocrtx':
    cfield = 'U'
    #cpal_fld = 'on2'
    cpal_fld = 'RdBu'
    tmin=-0.25 ;  tmax=0.25 ;  df = 0.05
    cunit = ''
    cb_jump = 1


else:
    print 'ERROR: we do not know cv!'
    sys.exit(0)



    


bt.chck4f(cf_lsm)
bt.chck4f(cf_in)
#id_fld = Dataset(cf_in)
#vtime = id_fld.variables['time_counter'][:]
#id_fld.close()
#Nt = len(vtime)

cv_lsm = 'tmask'

if cv_in == 'vozocrtx': cv_lsm = 'umask'


bt.chck4f(cf_lsm)
print '\n '+cv_lsm+' !!!'
id_lsm = Dataset(cf_lsm)
nb_dim = len(id_lsm.variables[cv_lsm].dimensions)
print ' The mesh_mask has '+str(nb_dim)+' dimmensions!'
if nb_dim==4: XMSK  = id_lsm.variables[cv_lsm][0,k1:k2,j1:j2,0]
if nb_dim==3: XMSK  = id_lsm.variables[cv_lsm][k1:k2,j1:j2,0]
if nb_dim==2: XMSK  = id_lsm.variables[cv_lsm][k1:k2,j1:j2]
(nk,nj) = nmp.shape(XMSK)
#XE1T2 = id_lsm.variables['e1t'][0,k1:k2,j1:j2]
#XE2T2 = id_lsm.variables['e2t'][0,k1:k2,j1:j2]
#(nk,nj) = nmp.shape(XE1T2)
#XE1T2 = XE1T2*XE1T2
#XE2T2 = XE2T2*XE2T2
id_lsm.close()

print 'Shape Arrays => nj,nk =', nj,nk

id_fld = Dataset(cf_in)
Nt = len(id_fld.variables[cv_in][:,0,0])
id_fld.close()

print ' *** Nt = ', Nt
print 'Done!\n'




pmsk = nmp.ma.masked_where(XMSK[:,:] > 0.2, XMSK[:,:]*0.+40.)



#font_rat
#params = { 'font.family':'Ubuntu',
params = { 'font.family':'Helvetica Neue',
           'font.weight':    'normal',
           'font.size':       int(9.*font_rat),
           'legend.fontsize': int(9.*font_rat),
           'xtick.labelsize': int(9.*font_rat),
           'ytick.labelsize': int(9.*font_rat),
           'axes.labelsize':  int(9.*font_rat) }
mpl.rcParams.update(params)
cfont_clb  = { 'fontname':'Helvetica Neue', 'fontweight':'medium', 'fontsize':int(8.*font_rat), 'color':'w'}
cfont_date = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(12.*font_rat), 'color':'w' }
cfont_mail = { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14.*font_rat), 'color':'0.8'}
cfont_titl = { 'fontname':'Helvetica Neue', 'fontweight':'light', 'fontsize':int(30.*font_rat), 'color':'w' }


# Colormaps for fields:
pal_fld = bcm.chose_colmap(cpal_fld)
if l_log_field:
    norm_fld = colors.LogNorm(  vmin = tmin, vmax = tmax, clip = False)
if l_pow_field:
    norm_fld = colors.PowerNorm(gamma=pow_field, vmin = tmin, vmax = tmax, clip = False)
else:
    norm_fld = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)


pal_lsm = bcm.chose_colmap('land_dark')
norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)




if cdt == '3h':
    dt = 3
elif cdt == '1h':
    dt = 1
else:
    print 'ERROR: unknown dt!'




print ' *** Dimension image:', rh*float(dpi), rh*yx_ratio*float(dpi),'\n'


ntpd = 24/dt


vm = vmn
if isleap(int(cyr0)): vm = vml
#print ' year is ', vm, nmp.sum(vm)

jd = int(cdd0) - 1
jm = int(cmn0)

for jt in range(jt0,Nt):

    jh = (jt*dt)%24
    jdc = (jt*dt)/24 + 1

    if jt%ntpd == 0: jd = jd + 1

    if jd == vm[jm-1]+1 and (jt)%ntpd == 0 :
        jd = 1
        jm = jm + 1

    ch = '%2.2i'%(jh)
    #cdc= '%3.3i'%(jdc)
    cd = '%3.3i'%(jd)
    cm = '%2.2i'%(jm)

    #print '\n\n *** jt, ch, cd, cm =>', jt, ch, cd, cm


    ct = str(datetime.datetime.strptime(cyr0+'-'+cm+'-'+cd+' '+ch, '%Y-%m-%j %H'))
    ct=ct[:5]+cm+ct[7:] #lolo bug !!! need to do that to get the month and not "01"
    print ' ct = ', ct
    cday  = ct[:10]   ; print ' *** cday  :', cday
    chour = ct[11:13] ; print ' *** chour :', chour



    cfig = 'figs/'+cv_in+'_NEMO_'+CNEMO+'_'+cbox+'_'+cday+'_'+chour+'_'+cpal_fld+'.'+fig_type

    fig = plt.figure(num = 1, figsize=(rh,rh*yx_ratio*3), dpi=None, facecolor='w', edgecolor='0.5')

    #ax  = plt.axes([0.065, 0.05, 0.9, 1.], axisbg = '0.5')
    ax  = plt.axes([0., 0., 1., 1.], axisbg = '0.5')

    vc_fld = nmp.arange(tmin, tmax + df, df)


    print "Reading record #"+str(jt)+" of "+cv_in+" in "+cf_in
    id_fld = Dataset(cf_in)
    XFLD  = id_fld.variables[cv_in][jt,k1:k2,j1:j2,0]
    id_fld.close()
    print "Done!"


    #if not l_show_lsm and jt == jt0: ( nk , nj ) = nmp.shape(XFLD)
    print '  *** dimension of array => ', nj, nk

    print "Ploting"
    cf = plt.imshow(XFLD[:,:], cmap = pal_fld, norm = norm_fld, interpolation='none')
    del XFLD
    print "Done!"

    #cm = plt.imshow(nmp.flipud(pmsk), cmap = pal_lsm, norm = norm_lsm, interpolation='none')

    #plt.axis([ 0, nj, 0, nk])

    #plt.title('NEMO: '+cfield+', coupled '+CNEMO+', '+cday+' '+chour+':00', **cfont_title)



    if l_show_cb:
        color_top='w'
        ax2 = plt.axes(vcb)
        clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend='both')
        if cb_jump > 1:
            cb_labs = [] ; cpt = 0
            for rr in vc_fld:
                if cpt % cb_jump == 0:
                    if df >= 1.: cb_labs.append(str(int(rr)))
                    if df <  1.: cb_labs.append(str(rr))
                else:
                    cb_labs.append(' ')
                    cpt = cpt + 1
                    clb.ax.set_xticklabels(cb_labs, **cfont_clb)
                    clb.set_label(cunit, **cfont_clb)
                    clb.ax.yaxis.set_tick_params(color=color_top) ; # set colorbar tick color
                    clb.outline.set_edgecolor(color_top) ; # set colorbar edgecolor
                    plt.setp(plt.getp(clb.ax.axes, 'xticklabels'), color=color_top) ; # set colorbar ticklabels

    del cf





    if l_show_date:
        xl = float(x_date)/rfact_zoom
        yl = float(y_date)/rfact_zoom
        ax.annotate('Date: '+cday+' '+chour+':00', xy=(1, 4), xytext=(xl,yl), **cfont_date)

    #ax.annotate('laurent.brodeau@ocean-next.fr', xy=(1, 4), xytext=(xl+150, 20), **cfont_mail)



    if l_annotate_name:
        xl = float(nxr)/20./rfact_zoom
        yl = float(nyr)/1.33/rfact_zoom
        ax.annotate(CNEMO, xy=(1, 4), xytext=(xl, yl), **cfont_titl)



    plt.savefig(cfig, dpi=dpi, orientation='portrait', facecolor='k')
    print cfig+' created!\n'
    plt.close(1)


    del cm, fig, ax
    if l_show_cb: del clb
