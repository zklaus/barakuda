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

import barakuda_colmap as bcm

import barakuda_tool as bt


#CNEMO = 'NATL60'
#CNEMO = 'NANUK025'

#color_top = 'white'
color_top = 'k'






    










fig_type='png'

narg = len(sys.argv)
if narg < 6: print 'Usage: '+sys.argv[0]+'<CONF> <file> <variable> <snapshot> <LSM_file>'; sys.exit(0)
CNEMO = sys.argv[1] ; cf_fld = sys.argv[2] ; cv_in=sys.argv[3] ; jt=int(sys.argv[4]) ; cf_lsm=sys.argv[5]


i2=0
j2=0

if CNEMO == 'NATL60':
    #i1 = 2880 ; j1 = 1060 ; rfact_zoom = 1. ; ny_res = 1080 ; nx_res = 1920 ; vcb = [0.01, 0.08, 0.98, 0.025]
    #i1 = 0 ; j1 = 150 ; rfact_zoom = 0.36  ; ny_res = 1080 ; nx_res = 1920 ; vcb = [0.01, 0.08, 0.98, 0.025]
    i1 = 4086 ; j1 = 603 ; i2 = 5421 ; j2 = 3119 ; rfact_zoom = 0.3319 ; vcb = [0.01, 0.08, 0.98, 0.025]

elif CNEMO == 'NANUK025':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.5, 0.875, 0.49, 0.02]

elif CNEMO == 'eNATL4':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.6, 0.11, 0.39, 0.025]

elif CNEMO == 'eNATL60':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.6, 0.1, 0.39, 0.025]

else:
    print '\n PROBLEM: "'+CNEMO+'" is an unknown config!!!'
    sys.exit(0)









l_log_field = False
cextend='both'

if cv_in in ['sosstsst','tos']:
    cfield = 'SST'
    tmin=0. ;  tmax=20.   ;  df = 1.
    cpal_fld = 'ncview_nrl'    
    cunit = r'SST ($^{\circ}$C)'
    cb_jump = 1


if cv_in in ['Bathymetry']:
    cfield = 'Bathymetry'
    tmin=100. ;  tmax=4500.   ;  df = 100.
    #cpal_fld = 'ocean'
    #cpal_fld = 'Blues'
    cpal_fld = 'PuBu'    
    cunit = r'Bathymetry (m)'
    cb_jump = 1
    l_log_field = True
    cextend='max'

    
if cv_in == 'sossheig':
    cfield = 'SSH'
    tmin=-0.3 ;  tmax=0.3   ;  df = 0.1
    cpal_fld = 'ncview_jaisnc'    
    cunit = r'SSH (m)'
    cb_jump = 1

elif cv_in == 'somxl010':
    cfield == 'MLD'
    tmin=50. ;  tmax=1500. ;  df = 50.
    cpal_fld = 'viridis_r'
    

# Time record stuff...
bt.chck4f(cf_fld)
id_fld = Dataset(cf_fld)
list_var = id_fld.variables.keys()
if 'time_counter' in list_var:    
    vtime = id_fld.variables['time_counter'][:]
    Nt = len(vtime)
    print '\n There is a "time_counter" in file '+cf_fld+' !'
    print '   => '+str(Nt)+' snapshots!'
    if jt <= 0 or jt > Nt: print ' PROBLEM: your time record does not exist!', jt ; sys.exit(0)
else:
    print '\n There is NO "time_counter" in file '+cf_fld+' !'
    Nt = 0    
id_fld.close()






bt.chck4f(cf_lsm)
print '\n *** Reading "tmask" in meshmask file...'
id_lsm = Dataset(cf_lsm)
Ni = id_lsm.dimensions['x'].size
Nj = id_lsm.dimensions['y'].size
if i2 == 0: i2 = Ni
if j2 == 0: j2 = Nj
XMSK  = id_lsm.variables['tmask'][0,0,j1:j2,i1:i2] ; # t, y, x
id_lsm.close()
print '      done.'


print '\n According to "tmask" the shape of the domain is Ni, Nj =', Ni, Nj


# Stuff for size of figure respecting pixels...
print '  *** we are going to show: i1,i2,j1,j2 =>', i1,i2,j1,j2, '\n'
nx_res = i2-i1
ny_res = j2-j1
yx_ratio = float(ny_res)/float(nx_res)
nxr = int(rfact_zoom*nx_res) ; # widt image (in pixels)
nyr = int(rfact_zoom*ny_res) ; # height image (in pixels)
dpi = 110
rh  = float(nxr)/float(dpi) ; # width of figure as for figure...
font_rat = nxr/1080.





pmsk = nmp.ma.masked_where(XMSK[:,:] > 0.2, XMSK[:,:]*0.+40.)





idx_oce = nmp.where(XMSK[:,:] > 0.5)

#font_rat
#params = { 'font.family':'Ubuntu',
params = { 'font.family':'Helvetica Neue',
           'font.weight':    'normal',
           'font.size':       int(12.*font_rat),
           'legend.fontsize': int(12.*font_rat),
           'xtick.labelsize': int(12.*font_rat),
           'ytick.labelsize': int(12.*font_rat),
           'axes.labelsize':  int(12.*font_rat) }
mpl.rcParams.update(params)
cfont_clb  = { 'fontname':'Helvetica Neue', 'fontweight':'medium', 'fontsize':int(12.*font_rat), 'color':color_top }
cfont_date = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(12.*font_rat), 'color':'w' }
cfont_mail = { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14.*font_rat), 'color':'0.8'}
cfont_titl = { 'fontname':'Helvetica Neue', 'fontweight':'light', 'fontsize':int(30.*font_rat), 'color':'w' }


# Colormaps for fields:
pal_fld = bcm.chose_colmap(cpal_fld)
if l_log_field:
    norm_fld = colors.LogNorm(  vmin = tmin, vmax = tmax, clip = False)
else:
    norm_fld = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)

pal_lsm = bcm.chose_colmap('land_dark')
norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)


if Nt == 0:
    cfig = cv_in+'_'+CNEMO+'_'+cpal_fld+'.'+fig_type    
else:
    cfig = 'snapshot_'+str(jt)+'_'+cv_in+'_'+CNEMO+'_'+cpal_fld+'.'+fig_type    

fig = plt.figure(num = 1, figsize=(rh,rh*yx_ratio), dpi=None, facecolor='w', edgecolor='0.5')

#ax  = plt.axes([0.065, 0.05, 0.9, 1.], axisbg = '0.5')
ax  = plt.axes([0., 0., 1., 1.], axisbg = '0.5')

vc_fld = nmp.arange(tmin, tmax + df, df)


print '\n *** Opening file '+cf_fld
id_fld = Dataset(cf_fld)
if Nt > 0:
    print '    => Reading record #'+str(jt)+' of '+cv_in+' in '+cf_fld
    XFLD  = id_fld.variables[cv_in][jt-1,j1:j2,i1:i2] ; # t, y, x
else:
    print '    => Reading 2D field '+cv_in+' in '+cf_fld+' (no time records...)'
    XFLD  = id_fld.variables[cv_in][j1:j2,i1:i2] ; # t, y, x

id_fld.close()
print '          Done!\n'


if XMSK.shape != XFLD.shape:
    print '\n PROBLEM: field and mask do not agree in shape!'
    print XMSK.shape , XFLD.shape
    sys.exit(0)

print '  *** Shape of field and mask => ', nmp.shape(XFLD)


del XMSK


print 'Ploting'
cf = plt.imshow(XFLD[:,:], cmap = pal_fld, norm = norm_fld, interpolation='none')
del XFLD
print 'Done!'



cm = plt.imshow(pmsk, cmap = pal_lsm, norm = norm_lsm, interpolation='none')

del pmsk


plt.axis([ 0, Ni, 0, Nj])

#plt.title('NEMO: '+cfield+', coupled '+CNEMO+', '+cday+' '+chour+':00', **cfont_title)


#ax2 = plt.axes([0.3, 0.08, 0.4, 0.025])

ax2 = plt.axes(vcb)
    
clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend=cextend)
if cb_jump > 1:
    cb_labs = [] ; cpt = 0
    for rr in vc_fld:
        if cpt % cb_jump == 0:
            if df >= 1.: cb_labs.append(str(int(rr)))
            if df <  1.: cb_labs.append(str(rr))
        else:
            cb_labs.append(' ')
        cpt = cpt + 1
    clb.ax.set_xticklabels(cb_labs)    
clb.set_label(cunit, **cfont_clb)
clb.ax.yaxis.set_tick_params(color=color_top) ; # set colorbar tick color    
clb.outline.set_edgecolor(color_top) ; # set colorbar edgecolor         
plt.setp(plt.getp(clb.ax.axes, 'xticklabels'), color=color_top) ; # set colorbar ticklabels
    
del cf


#x_annot = nxr-nxr*0.22.*font_rat ; y_annot = 150
x_annot = 650 ; y_annot = 1035

#ax.annotate('Date: '+cday+' '+chour+':00',   xy=(1, 4), xytext=(x_annot,    y_annot), **cfont_date)

#ax.annotate('laurent.brodeau@ocean-next.fr', xy=(1, 4), xytext=(x_annot+150, 20), **cfont_mail)


xl = float(nxr)/20./rfact_zoom
yl = float(nyr)/1.2/rfact_zoom
ax.annotate(CNEMO, xy=(1, 4), xytext=(xl, yl), **cfont_titl)



plt.savefig(cfig, dpi=dpi, orientation='portrait', facecolor='k')
print cfig+' created!\n'
plt.close(1)


del cm, fig, ax, clb


