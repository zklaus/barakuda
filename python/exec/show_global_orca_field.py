#!/usr/bin/env python

#    B a r a K u d a
#
#
#    L. Brodeau, 2017

import sys
import os
from string import replace
import numpy as nmp

from netCDF4 import Dataset

from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.font_manager as font_manager

import warnings
warnings.filterwarnings("ignore")

import barakuda_orca   as bo
import barakuda_ncio   as bnc
import barakuda_colmap as bcm

import barakuda_tool as bt


l_force_mask = True
nx_exclude_north = 300 ; # nb points to exclude at the north...
nx_exclude_south = 200 ; #               ''            south...

lon_reorg = True
l_show_colorbar = False

color_text_colorbar = 'w'
color_top = 'k'
color_continents='#9C5536'


#nx_res = 1600
#nx_res = 1920 ; # number of pixels you want in the final image...
nx_res = 4680 ; # full res in our case!!!
rDPI=100.
rh = float(nx_res)/rDPI

longitude = 50. ; #
latitude  = 0. ; # fixed latitude to view from

#lforce_mask = True
#year_ref_ini = 1990

cf_mm = '/data/gcm_output/NEMO/ORCA12.L75/ORCA12.L75-I/mesh_mask.nc4' ; # NEMO mesh_mask file...

fig_type='png'

narg = len(sys.argv)
if narg < 4: print 'Usage: '+sys.argv[0]+' <file> <variable> <# snapshot>'; sys.exit(0)
cf_in = sys.argv[1] ; cv_in=sys.argv[2] ; jt=int(sys.argv[3])-1





if cv_in == 'tap':
    cfield = 'TAP'
    tmin=0. ;  tmax=600.   ;  dtemp = 10.    
    #cpal_fld = 'ncview_hotres'
    #cpal_fld = 'hot_r'
    cpal_fld = 'tap' ; log_ctrl=0.1  ; # 0 if you want no log involved in the colorscale...    
    cunit = '(MW)'
    cb_jump = 4
elif cv_in == 'sosstsst':
    cfield = 'SST'
    log_ctrl=0  ; # 0 if you want no log involved in the colorscale...
    tmin=-2. ;  tmax=32.   ;  dtemp = 2.
    #cpal_fld = 'ncview_hotres'
    #cpal_fld = 'hot_r'
    cpal_fld = 'tap'
    cunit = r'deg.C'
    cb_jump = 2




j1=nx_exclude_north

    
bt.chck4f(cf_in)

bt.chck4f(cf_mm)
id_mm = Dataset(cf_mm)
XLONo  =  id_mm.variables['glamt'][0,:,:]
(njo,nio) = nmp.shape(XLONo)
j1=nx_exclude_south
j2=njo-nx_exclude_north
del XLONo
XLONo  =  id_mm.variables['glamt'][0,j1:j2,:]
XLATo  =  id_mm.variables['gphit'][0,j1:j2,:]
if l_force_mask: XMSKo  = id_mm.variables['tmask'][0,0,j1:j2,:]
id_mm.close()




#idx_oce = nmp.where(XMSK[:,:] > 0.5)

font_corr = float(nx_res)/1200.


params = { 'font.family':'Helvetica Neue',
           'font.weight':    'light',
           'font.size':       int(12*font_corr),
           'legend.fontsize': int(12*font_corr),
           'xtick.labelsize': int(12*font_corr),
           'ytick.labelsize': int(12*font_corr),
           'axes.labelsize':  int(12*font_corr) }
mpl.rcParams.update(params)


#path = '/home/laurent/.fonts/OSX_conv_LinuX/HelveticaNeue.ttf'
#prop = font_manager.FontProperties(fname=path)
#prop.set_weight = 'light'
#mpl.rcParams['font.family'] = prop.get_name()
#mpl.rcParams['font.weight'] = 'light'

#cfont_clb   = { 'fontweight':'ultra-light', 'fontsize':int(13*font_corr), 'color':'white' }
#cfont_title = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(18*font_corr), 'color':'white' }
#cfont_mail  = { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(10*font_corr), 'color':'0.5'}





# Colormaps for fields:
pal_fld = bcm.chose_colmap(cpal_fld, log_ctrl=log_ctrl)
norm_fld = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)

vc_fld = nmp.arange(tmin, tmax + dtemp, dtemp)


#pal_ice = bcm.chose_colmap(cpal_ice)
#norm_ice = colors.Normalize(vmin = rmin_ice, vmax = 1, clip = False)

#pal_mm = bcm.chose_colmap('blk')
#norm_mm = colors.Normalize(vmin = 0., vmax = 1., clip = False)



print ''



vproj = [ 'ortho' ]
#vproj = [ 'gall', 'mill'  ,'merc'  ,'cyl'  ,'mbtfpq'  ,'kav7'  ,'moll'  ,'eck4'  ,'robin'  ,'hammer' ]
#vproj = [ 'mbtfpq' ]

#cd = str(datetime.datetime.strptime('1990 '+ct, '%Y %j'))
#cdate = cd[:10] ; print ' *** Date :', cdate

id_fld = Dataset(cf_in)
XFLDo = id_fld.variables[cv_in][jt,j1:j2,:]
id_fld.close()
print "  => done!"


if lon_reorg:
    XLAT = bo.lon_reorg_orca(XLATo, XLONo, ilon_ext=30)
    XFLD = bo.lon_reorg_orca(XFLDo, XLONo, ilon_ext=30)
    XLON = bo.lon_reorg_orca(XLONo, XLONo, ilon_ext=30)
    if l_force_mask: XMSK = bo.lon_reorg_orca(XMSKo, XLONo, ilon_ext=30)
    print "shape old XLON =>", nmp.shape(XLONo)
    print "shape new XLON =>", nmp.shape(XLON)
    #idn = nmp.where( XLONa < 0. )
    #XLONa[idn] = XLONa[idn] + 360.
    #bnc.write_2d_mask('zshowa.nc', XFLDa, xlon=XLONa, xlat=XLATa, name='tap')
    #XLAT = bo.lon_reorg_orca(XLATa, XLONa, v_junc_i_p=355., v_junc_i_m=5.)
    #XFLD = bo.lon_reorg_orca(XFLDa, XLONa, v_junc_i_p=355., v_junc_i_m=5.)
    #XLON = bo.lon_reorg_orca(XLONa, XLONa, v_junc_i_p=355., v_junc_i_m=5.)
    #print "shape new XLON =>", nmp.shape(XLON)
    #del XLONa, XLATa, XFLDa
    
else:
    XLAT = XLATo
    XFLD = XFLDo
    XLON = XLONo
    XMSK = XMSKo

del XLATo,XLONo,XFLDo,XMSKo

(nj,ni) = nmp.shape(XLON)

#lolo: bnc.write_2d_mask('zshow.nc', XFLD, xlon=XLON, xlat=XLAT, name='tap')

pmsk = nmp.ma.masked_where(XMSK[:,:] > 0.2, XMSK[:,:]*0.+40.)



# No proj :
#pal_lsm = bcm.chose_colmap('terre')
pal_lsm = bcm.chose_colmap('land')
norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)



# Creating colorbar in a dfferent image:
fig = plt.figure(num = 2, figsize=(rh,rh/18.), dpi=None) #, facecolor='w', edgecolor='0.5')
ax2 = plt.axes([0., 0., 1., 1.], axisbg = None)
ax2 = plt.axes([0.2, 0.5, 0.6, 0.4])
clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend='both')
cb_labs = [] ; cpt = 0
for rr in vc_fld:
    if cpt % cb_jump == 0:
        cb_labs.append(str(int(rr)))
    else:
        cb_labs.append(' ')
    cpt = cpt + 1
clb.ax.set_xticklabels(cb_labs)
clb.set_label(cunit, color=color_text_colorbar)
clb.ax.yaxis.set_tick_params(color=color_top) ; # set colorbar tick color    
clb.outline.set_edgecolor(color_top) ; # set colorbar edgecolor         
plt.setp(plt.getp(clb.ax.axes, 'xticklabels'), color=color_top) ; # set colorbar ticklabels

cbytick_obj = plt.getp(clb.ax.axes, 'xticklabels')                #tricky
plt.setp(cbytick_obj, color=color_text_colorbar)

#for l in clb.ax.xaxis.get_ticklabels():
#    l.set_family("Fixed")


plt.savefig('colorbar.svg', dpi=rDPI, orientation='portrait', transparent=True)




fig = plt.figure(num = 1, figsize=(rh,rh*float(nj)/float(ni)), dpi=None) #, facecolor='w', edgecolor='0.5')

#ax  = plt.axes([0.065, 0.05, 0.9, 1.], axisbg = '0.5')
ax  = plt.axes([0., 0., 1., 1.], axisbg = '0.5')








cf = plt.imshow(XFLD[:,:] , cmap = pal_fld, norm = norm_fld)

cm = plt.imshow(pmsk, cmap = pal_lsm, norm = norm_lsm)
plt.axis([ 0, ni, 0, nj])


if l_show_colorbar:
    ax2 = plt.axes([0.2, 0.08, 0.6, 0.025])
    clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend='both')
    cb_labs = [] ; cpt = 0
    for rr in vc_fld:
        if cpt % cb_jump == 0:
            cb_labs.append(str(int(rr)))
        else:
            cb_labs.append(' ')
        cpt = cpt + 1
    clb.ax.set_xticklabels(cb_labs)
    clb.set_label(cunit, **cfont_clb)
    clb.ax.yaxis.set_tick_params(color=color_top) ; # set colorbar tick color    
    clb.outline.set_edgecolor(color_top) ; # set colorbar edgecolor         
    plt.setp(plt.getp(clb.ax.axes, 'xticklabels'), color=color_top) ; # set colorbar ticklabels
        
del cf
    

#ax.annotate('Date: '+cday+' '+chour+':00', xy=(1, 4), xytext=(nx_res-nx_res*0.22, 95), **cfont_title)

#ax.annotate('laurent.brodeau@ocean-next.fr', xy=(1, 4), xytext=(nx_res-nx_res*0.2, 20), **cfont_mail)

plt.savefig('figure01.png', dpi=rDPI, orientation='portrait', transparent=True) ; #facecolor='k')






#lulu





sys.exit(0)

    
for cproj in vproj:

    cfig = 'figs/'+cv_in+'_NEMO'+'_'+cproj+'.'+fig_type

    print ' *** reference longitude =', longitude

    #fig = plt.figure(num = 1, figsize=(rh,1.167*rh), dpi=None, facecolor='b', edgecolor='k')
    #ax  = plt.axes([0.005, 0.05, 0.99, 0.99], axisbg = '0.35')
    fig = plt.figure(num = 1, figsize=(rh,0.7*rh), dpi=None, facecolor='b', edgecolor='k')
    ax  = plt.axes([0.05, 0.05, 0.99, 0.99], axisbg = '0.35')

    plt.title('ORCA12: '+cfield, **cfont_title)

    print ' *** Creating projection '+cproj
    if cproj == 'ortho':
        carte = Basemap(projection=cproj, lat_0=latitude, lon_0=longitude, resolution='h')
    else:
        carte = Basemap(projection=cproj, lon_0=longitude, resolution='c')

    x0,y0 = carte(XLON,XLAT)

    #if l_force_mask: XFLD = nmp.ma.masked_where(XMSK[:,:] < 0.5, XFLD[:,:])

    print ' *** Ploting on map...'
    cf = carte.pcolor(x0, y0, XFLD, cmap=pal_fld, norm=norm_fld)

    carte.drawcoastlines()
    carte.fillcontinents(color='grey')
    #carte.drawmapboundary()
    #carte.drawmeridians(nmp.arange(-180,180,20), labels=[0,0,0,1])
    #carte.drawparallels(nmp.arange( -90, 90,10), labels=[1,0,0,0])

    plt.annotate('L. Brodeau, brodeau@gmail.com', xy=(0.7, 0.1), xycoords='figure fraction', **cfont_mail)

    # Colorbar:
    ax2 = plt.axes([0.005, 0.06, 0.99, 0.025])
    clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend='both')
    if cb_jump > 1:
        cb_labs = [] ; cpt = 0
        for rr in vc_fld:
            if cpt % cb_jump == 0:
                cb_labs.append(str(int(rr)))
            else:
                cb_labs.append(' ')
            cpt = cpt + 1
    clb.ax.set_xticklabels(cb_labs)
    clb.set_label(cunit, **cfont_clb)
    clb.ax.yaxis.set_tick_params(color='white') ; # set colorbar tick color
    clb.outline.set_edgecolor('white') ; # set colorbar edgecolor
    plt.setp(plt.getp(clb.ax.axes, 'xticklabels'), color='white') ; # set colorbar ticklabels

    print ' *** Saving figure...'
    #plt.savefig(cfig, dpi=160, orientation='portrait', facecolor='#06051F')
    plt.savefig(cfig, dpi=100, orientation='portrait', facecolor='k')
    print '  => '+cfig+' created!\n'
    plt.close(1)

