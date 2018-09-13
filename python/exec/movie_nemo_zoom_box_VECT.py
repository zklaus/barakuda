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

#CNEMO = 'eNATL60'
#CNEMO = 'NATL60'
#CNEMO = 'NANUK025'

color_top = 'white'
#color_top = 'k'



cv_out = 'unknown'

#jt0 = 248
jt0 = 0


i2=0
j2=0
l_show_lsm = True
l_do_ice  = True
l_show_cb = True
l_show_date = True
l_log_field = False
l_pow_field = False
l_annotate_name = True

l_do_curl = True
romega = 2.*nmp.pi/86400.0





narg = len(sys.argv)
if narg < 8: print 'Usage: '+sys.argv[0]+' <NEMOCONF> <fileX> <varX> <fileY> <varY> <LSM_file> <YYYYMMDD (start)>'; sys.exit(0)
CNEMO  = sys.argv[1]
cfx_in = sys.argv[2] ; cvx_in = sys.argv[3]
cfy_in = sys.argv[4] ; cvy_in = sys.argv[5]
cf_lsm = sys.argv[6] ; cf_date0=sys.argv[7]







if CNEMO == 'eNATL60':
    Ni0 = 8354-1
    Nj0 = 4729-1
    l_do_ice  = False
    l_show_cb = False
    l_show_date = True
    cdt= '1h'
    cbox = 'FullMed' ; i1=5400 ; j1=1530 ; i2=Ni0 ; j2=3310 ; rfact_zoom = 0.79 ; vcb = [0.5, 0.875, 0.485, 0.02] ; font_rat = 2.*rfact_zoom ; l_annotate_name=False
    #cbox = 'ALL' ; i1=0 ; j1=0 ; i2=Ni0 ; j2=Nj0 ; rfact_zoom = 0.3047 ; vcb = [0.59, 0.1, 0.38, 0.018] ; font_rat = 8.*rfact_zoom
    x_date = 1900 ; y_date = 20 ; # where to put the date

if CNEMO == 'NATL60':
    Ni0 = 5422-1
    Nj0 = 3454-1
    #l_pow_field = True ; pow_field = 1.5
    l_do_ice  = False
    l_show_cb = False
    l_show_date = False
    cdt = '1h'
    #cbox = 'zoom1' ; i1 = 1800 ; j1 = 950 ; i2 = i1+1920 ; j2 = j1+1080 ; rfact_zoom = 1. ; vcb = [0.5, 0.875, 0.485, 0.02] ; font_rat = 8.*rfact_zoom ; l_show_lsm = False
    #cbox = 'zoom1' ; i1 = 1800 ; j1 = 950 ; i2 = i1+2560 ; j2 = j1+1440 ; rfact_zoom = 1. ; vcb = [0.5, 0.875, 0.485, 0.02] ; font_rat = 8.*rfact_zoom
    cbox = 'ALL' ; i1=0 ; j1=0 ; i2=Ni0 ; j2=Nj0 ; rfact_zoom = 0.4 ; vcb = [0.59, 0.1, 0.38, 0.018] ; font_rat = 4.*rfact_zoom
    x_date = 350 ; y_date = 7 ; # where to put the date


if CNEMO == 'NANUK025':
    l_do_ice = True
    cdt = '3h'; cbox = 'ALL' ; i1 = 0 ; j1 = 0 ; i2 = 492 ; j2 = 614 ; rfact_zoom = 2. ; vcb = [0.5, 0.875, 0.485, 0.02] ; font_rat = 8.*rfact_zoom
    x_date = 350 ; y_date = 7 ; # where to put the date


print '\n rfact_zoom = ', rfact_zoom
print ' font_rat = ', font_rat, '\n'

nx_res = i2-i1
ny_res = j2-j1

print ' *** nx_res, ny_res =', nx_res, ny_res


print ' i1,i2,j1,j2 =>', i1,i2,j1,j2

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

# Ice:

if l_do_curl: cv_out = 'RV'


if l_do_ice:
    cv_ice  = 'siconc'
    cf_ice = replace(cfx_in, 'grid_T', 'icemod')
    rmin_ice = 0.5
    cpal_ice = 'ncview_bw'
    vcont_ice = nmp.arange(rmin_ice, 1.05, 0.05)

if cvx_in in ['sosstsst','tos']:
    cfield = 'SST'
    tmin=0. ;  tmax=30.   ;  df = 2. ; cpal_fld = 'ncview_nrl'
    #tmin=0. ;  tmax=32.   ;  df = 2. ; cpal_fld = 'viridis'
    #tmin=4. ;  tmax=20.   ;  df = 1. ; cpal_fld = 'PuBu'
    cunit = r'SST ($^{\circ}$C)'
    cb_jump = 2
    l_show_cb = True

if cvx_in == 'sossheig':
    cfield = 'SSH'
    #tmin=-0.5 ;  tmax=0.5   ;  df = 0.05
    #tmin=-1.2 ;  tmax=2.3   ;  df = 0.05
    #cpal_fld = 'ncview_jaisnc'
    #cpal_fld = 'PuBu'
    #cpal_fld = 'RdBu'
    #cpal_fld = 'BrBG'
    #
    #cpal_fld = 'on3' ; tmin=-1.2 ;  tmax=2.3   ;  df = 0.05
    cpal_fld = 'on2' ; tmin=-1.2 ;  tmax=1.2   ;  df = 0.05
    #cpal_fld = 'coolwarm' ; tmin=-1. ;  tmax=1.   ;  df = 0.05 
    #cpal_fld = 'RdBu_r' ; tmin=-0.9 ;  tmax=-tmin   ;  df = 0.05 
    #cpal_fld = 'gray_r' ; tmin=-0.3 ;  tmax=0.3   ;  df = 0.05
    #cpal_fld = 'bone_r' ; tmin=-0.9 ;  tmax=-tmin   ;  df = 0.05 
    cunit = r'SSH (m)'
    cb_jump = 1

elif cvx_in=='sozocrtx' and cvy_in=='somecrty':
    cfield = 'RV'
    cpal_fld = 'on2' ; tmin=-1. ;  tmax=1. ;  df = 0.05
    cunit = ''
    cb_jump = 1

elif cvx_in=='vozocrtx' and cvy_in=='vomecrty':
    #l_3d_field = True
    cfield = 'RV'
    cpal_fld = 'on2' ; tmin=-1. ;  tmax=1. ;  df = 0.05
    cunit = ''
    cb_jump = 1


else:
    print 'ERROR: we do not know cvx_in and cvy_in!'
    sys.exit(0)



    

if l_do_ice: bt.chck4f(cf_ice)

bt.chck4f(cf_lsm)
bt.chck4f(cfx_in)
bt.chck4f(cfy_in)

id_fx = Dataset(cfx_in)
vtime = id_fx.variables['time_counter'][:]
id_fx.close()

Nt = len(vtime)

if l_show_lsm or l_do_curl:
    print "\nReading record metrics in "+cf_lsm
    id_lsm = Dataset(cf_lsm)
    nb_dim = len(id_lsm.variables['tmask'].dimensions)
    print ' The mesh_mask has '+str(nb_dim)+' dimmensions!'
    if l_show_lsm:
        if nb_dim==4: XMSK = id_lsm.variables['tmask'][0,0,j1:j2,i1:i2]
        if nb_dim==3: XMSK = id_lsm.variables['tmask'][0,j1:j2,i1:i2]
        if nb_dim==2: XMSK = id_lsm.variables['tmask'][j1:j2,i1:i2]
    if l_do_curl:
        # e2v, e1u, e1f, e2f
        e2v = id_lsm.variables['e2v'][0,j1:j2,i1:i2]
        e1u = id_lsm.variables['e1u'][0,j1:j2,i1:i2]
        e1f = id_lsm.variables['e1f'][0,j1:j2,i1:i2]
        e2f = id_lsm.variables['e2f'][0,j1:j2,i1:i2]
        ff  = id_lsm.variables['gphif'][0,j1:j2,i1:i2]
        #ff  = id_lsm.variables['ff'][0,j1:j2,i1:i2]
        if nb_dim==4: XMSK = id_lsm.variables['fmask'][0,0,j1:j2,i1:i2]
        if nb_dim==3: XMSK = id_lsm.variables['fmask'][0,j1:j2,i1:i2]
        if nb_dim==2: XMSK = id_lsm.variables['fmask'][j1:j2,i1:i2]
        # Coriolis Parameter:
        ff[:,:] = 2.*romega*nmp.sin(ff[:,:]*nmp.pi/180.0)        
    (nj,ni) = nmp.shape(XMSK)
    id_lsm.close()

    print 'Shape Arrays => ni,nj =', ni,nj
    
    print 'Done!\n'

    
if l_show_lsm: pmsk = nmp.ma.masked_where(XMSK[:,:] > 0.2, XMSK[:,:]*0.+40.)



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


if l_show_lsm:
    pal_lsm = bcm.chose_colmap('land_dark')
    norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)

if l_do_ice:
    pal_ice = bcm.chose_colmap(cpal_ice)
    norm_ice = colors.Normalize(vmin = rmin_ice, vmax = 1, clip = False)



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


Xplot = nmp.zeros((nj,ni))

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



    cfig = 'figs/'+cv_out+'_NEMO_'+CNEMO+'_'+cbox+'_'+cday+'_'+chour+'_'+cpal_fld+'.'+fig_type

    fig = plt.figure(num = 1, figsize=(rh,rh*yx_ratio), dpi=None, facecolor='w', edgecolor='0.5')

    ax  = plt.axes([0., 0., 1., 1.], axisbg = '0.5')

    vc_fld = nmp.arange(tmin, tmax + df, df)


    print "Reading record #"+str(jt)+" of "+cvx_in+" in "+cfx_in
    id_fx = Dataset(cfx_in)
    if not l_3d_field:
        XFLD  = id_fx.variables[cvx_in][jt,j1:j2,i1:i2] ; # t, y, x
    else:
        print 'j1:j2 =', j1,j2
        print 'i1:i2 =', i1,i2
        XFLD  = id_fx.variables[cvx_in][jt,0,j1:j2,i1:i2] ; # t, y, x
    id_fx.close()
    print "Done!"
    print "Reading record #"+str(jt)+" of "+cvy_in+" in "+cfy_in
    id_fy = Dataset(cfy_in)
    if not l_3d_field:
        YFLD  = id_fy.variables[cvy_in][jt,j1:j2,i1:i2] ; # t, y, x
    else:
        YFLD  = id_fy.variables[cvy_in][jt,0,j1:j2,i1:i2] ; # t, y, x
    id_fy.close()
    print "Done!"

    
    if l_do_curl:
        
        print '\nComputing curl...'
        lx = nmp.zeros((nj,ni))
        ly = nmp.zeros((nj,ni))
        
        lx[:,1:ni-1] =   e2v[:,2:ni]*YFLD[:,2:ni] - e2v[:,1:ni-1]*YFLD[:,1:ni-1] 
        ly[1:nj-1,:] = - e1u[2:nj,:]*XFLD[2:nj,:] + e1u[1:nj-1,:]*XFLD[1:nj-1,:]

        Xplot[:,:] = ( lx[:,:] + ly[:,:] )*XMSK[:,:] / ( e1f[:,:]*e2f[:,:]*ff[:,:] )         # Relative Vorticity...
        
        ##XFLD[:,:] = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ) - e2v(ji,jj) * vn(ji,jj)    &
        ##                                      &          - e1u(ji  ,jj+1) * un(ji  ,jj+1) + e1u(ji,jj) * un(ji,jj)  ) &
        ##                    &          * fmask(ji,jj) / ( e1f(ji,jj) * e2f(ji,jj) )

        del lx, ly
        print '... curl computed!\n'
        

    #Xplot[:,:] = YFLD[:,:]
    #Xplot[:,:] = XFLD[:,:]
    #Xplot[:,:] = ff[:,:]
    #if not l_show_lsm and jt == jt0: ( nj , ni ) = nmp.shape(XFLD)
    #print '  *** dimension of array => ', ni, nj


    del XFLD,YFLD
    
    print "Ploting"
    
    cf = plt.imshow(Xplot[:,:], cmap = pal_fld, norm = norm_fld, interpolation='none')

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
        ci = plt.imshow(pice, cmap = pal_ice, norm = norm_ice, interpolation='none') ; del pice, ci
        del XICE


    if l_show_lsm: cm = plt.imshow(pmsk, cmap = pal_lsm, norm = norm_lsm, interpolation='none')

    plt.axis([ 0, ni, 0, nj])

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
