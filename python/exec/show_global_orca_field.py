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
import barakuda_plot   as bp

import barakuda_tool as bt

CONF = 'ORCA12'


ldebug = False

l_add_lon_lat = False
l_draw_lat_lon_only = False
l_force_mask = True

lon_reorg = True
l_show_colorbar = False

ny_exclude_north = 0 ; # nb points to exclude at the north...
ny_exclude_south = 0 ; #               ''            south...
nx_exclude_west  = 0 ; #               ''            west



if CONF == 'ORCA12':
    ny_exclude_north = 350 ; # nb points to exclude at the north...
    ny_exclude_south = 150 ; #               ''            south...
    nx_exclude_west  = 300 ; #               ''            west
    cpal_fld = 'tap1'
    #cpal_fld = 'tap2'



    

color_text_colorbar = 'k'
color_stff_colorbar = 'k'
#color_continents    = '#9C5536'
#color_continents    = '#EDD0AB'
color_continents    = '0.75'




log_ctrl=0.1  ; # 0 if you want no log involved in the colorscale...
#b: log_ctrl=0.05  ; # 0 if you want no log involved in the colorscale...    


ilon_ext = 35 ; # Map extension on the RHS in degrees....

#nx_res = 1600
#nx_res = 1920 ; # number of pixels you want in the final image...
nx_res = -1 ; # trigers full res!!!
rDPI=100.


longitude = 50. ; #
latitude  = 0. ; # fixed latitude to view from

#lforce_mask = True
#year_ref_ini = 1990

cf_mm = '/data/gcm_output/NEMO/ORCA12.L75/ORCA12.L75-I/mesh_mask.nc4' ; # NEMO mesh_mask file...

fig_type='png'
if l_draw_lat_lon_only: fig_type='svg'

narg = len(sys.argv)
if narg < 4: print 'Usage: '+sys.argv[0]+' <file> <variable> <# snapshot>'; sys.exit(0)
cf_in = sys.argv[1] ; cv_in=sys.argv[2] ; jt=int(sys.argv[3])-1





if cv_in == 'tap':
    cfield = 'TAP'
    tmin=0. ;  tmax=600.   ;  dtemp = 10.    
    #cpal_fld = 'ncview_hotres'
    #cpal_fld = 'hot_r'
    #cpal_fld = 'tap1' ; log_ctrl=0.1  ; # 0 if you want no log involved in the colorscale...    
    cunit = '(MW)'
    cb_jump = 4
elif cv_in == 'sosstsst':
    cfield = 'SST'
    log_ctrl=0  ; # 0 if you want no log involved in the colorscale...
    tmin=-2. ;  tmax=32.   ;  dtemp = 2.
    #cpal_fld = 'ncview_hotres'
    #cpal_fld = 'hot_r'
    cunit = r'deg.C'
    cb_jump = 2




j1=ny_exclude_north

    
bt.chck4f(cf_in)

bt.chck4f(cf_mm)
id_mm = Dataset(cf_mm)
XLONo  =  id_mm.variables['glamt'][0,:,:]
(njo,nio) = nmp.shape(XLONo)
j1=ny_exclude_south
j2=njo-ny_exclude_north
del XLONo
XLONo  =  id_mm.variables['glamt'][0,j1:j2,:]
XLATo  =  id_mm.variables['gphit'][0,j1:j2,:]
if l_force_mask: XMSKo  = id_mm.variables['tmask'][0,0,j1:j2,:]
id_mm.close()







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
    XLAT = bo.lon_reorg_orca(XLATo, XLONo, ilon_ext=ilon_ext)
    XFLD = bo.lon_reorg_orca(XFLDo, XLONo, ilon_ext=ilon_ext)
    XLON = bo.lon_reorg_orca(XLONo, XLONo, ilon_ext=ilon_ext)
    if l_force_mask: XMSK = bo.lon_reorg_orca(XMSKo, XLONo, ilon_ext=ilon_ext)
    print "shape original XLON =>", nmp.shape(XLONo)
    print "shape intermediate XLON =>", nmp.shape(XLON)    
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


# Do we remove some RHS points:
XLON = XLON[:,nx_exclude_west:]
XLAT = XLAT[:,nx_exclude_west:]
XMSK = XMSK[:,nx_exclude_west:]
XFLD = XFLD[:,nx_exclude_west:]


(nj,ni) = nmp.shape(XLON)


print "shape final XLON =>", (nj,ni)


if nx_res == -1: nx_res = ni

#lolo: bnc.write_2d_mask('zshow.nc', XFLD, xlon=XLON, xlat=XLAT, name='tap')








# Time for figure...

font_corr = float(nx_res)/1200.

rh = float(nx_res)/rDPI

params = { 'font.family':'Helvetica Neue',
           'font.weight':    'light',
           'font.size':       int(14*font_corr),
           'legend.fontsize': int(14*font_corr),
           'xtick.labelsize': int(14*font_corr),
           'ytick.labelsize': int(14*font_corr),
           'axes.labelsize':  int(14*font_corr) }
mpl.rcParams.update(params)





# Creating colorbar in a dfferent image:
fig = plt.figure(num = 2, figsize=(rh,rh/18.), dpi=None) #, facecolor='w', edgecolor='0.5')
ax2 = plt.axes([0., 0., 1., 1.], facecolor = None)
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
clb.ax.yaxis.set_tick_params(color=color_stff_colorbar) ; # set colorbar tick color    
clb.outline.set_edgecolor(color_stff_colorbar) ; # set colorbar edgecolor         
plt.setp(plt.getp(clb.ax.axes, 'xticklabels'), color=color_stff_colorbar) ; # set colorbar ticklabels

cbytick_obj = plt.getp(clb.ax.axes, 'xticklabels')                #tricky
plt.setp(cbytick_obj, color=color_text_colorbar)

#for l in clb.ax.xaxis.get_ticklabels():
#    l.set_family("Fixed")


plt.savefig('colorbar_p'+cpal_fld+'_cc'+color_continents+'.svg', dpi=rDPI, orientation='portrait', transparent=True)




##################################################


fig = plt.figure(num = 1, figsize=(rh,rh*float(nj)/float(ni)), dpi=None) #, facecolor='r', edgecolor='0.5') #, facecolor='w', edgecolor='0.5')

ax  = plt.axes([0., 0., 1., 1.], facecolor = color_continents)
if ldebug: ax  = plt.axes([0.1, 0.1, 0.8, 0.8], facecolor = 'w')



if l_add_lon_lat or l_draw_lat_lon_only:
    VX = nmp.zeros(ni) ; VY = nmp.zeros(nj)
    VX[:] = XLON[100,:]
    dlong  = abs(VX[11] - VX[10])
    VX0 = nmp.arange(nx_exclude_west,ni+nx_exclude_west)  # lulu????
    VX = VX0*dlong + dlong/2.
    print VX,'\n'

    VY[:] = XLAT[:,nmp.argmax(XLAT[nj-1,:])]
    print VY,'\n'

    if not l_draw_lat_lon_only:
        cf = plt.pcolormesh(VX, VY, XFLD, cmap = pal_fld, norm = norm_fld)
        # Masking:
        pal_lsm = bcm.chose_colmap('terre')
        pal_lsm = bcm.chose_colmap('land')
        norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)
        pmsk = nmp.ma.masked_where(XMSK[:,:] > 0.2, XMSK[:,:]*0.+40.)
        cm = plt.pcolormesh(VX, VY, pmsk, cmap = pal_lsm, norm = norm_lsm)

    
    [vvx, vvy, clon, clat] = bp.__name_coor_ticks__(lon_ext=ilon_ext);
    plt.yticks(vvy,clat) ; plt.xticks(vvx,clon)
    plt.axis([ min(VX), 360.+ilon_ext-2., min(VY), max(VY)])

    
else:
    # Masking:
    idx_land = nmp.where(XMSK[:,:] < 0.2)
    XFLD[idx_land] = nmp.nan
    cf = plt.imshow(XFLD[:,:] , cmap = pal_fld, norm = norm_fld)
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
    clb.ax.yaxis.set_tick_params(color=color_stff_colorbar) ; # set colorbar tick color    
    clb.outline.set_edgecolor(color_stff_colorbar) ; # set colorbar edgecolor         
    plt.setp(plt.getp(clb.ax.axes, 'xticklabels'), color=color_stff_colorbar) ; # set colorbar ticklabels
        
#del cf   

#ax.annotate('Date: '+cday+' '+chour+':00', xy=(1, 4), xytext=(nx_res-nx_res*0.22, 95), **cfont_title)

#ax.annotate('laurent.brodeau@ocean-next.fr', xy=(1, 4), xytext=(nx_res-nx_res*0.2, 20), **cfont_mail)

plt.savefig('figure01_p'+cpal_fld+'_cc'+color_continents+'.'+fig_type, dpi=rDPI, orientation='portrait') #, transparent=True) ; #facecolor='k')




