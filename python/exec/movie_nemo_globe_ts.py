#!/usr/bin/env python

#       B a r a K u d a
#
#
#    L. Brodeau, 2017

import sys
import os
from string import replace
import numpy as nmp

from netCDF4 import Dataset

from mpl_toolkits.basemap import Basemap
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

#nsts = 4 ; # sub time-steping for creating interpolated frames (smoother video)
nsts = 1 ; # sub time-steping for creating interpolated frames (smoother video)

long_start = 0 ; # longitude to start the movie from...

latitude = 15. ; # fixed latitude to view from

lforce_mask = True
year_ref_ini = 1990

jt0 = 0
jtN = 0
#jt0 = 0

#CTATM = 'T255'
CTATM = 'T1279'

cbox = 'Pac'

fig_type='png'

narg = len(sys.argv)
if narg < 6: print 'Usage: '+sys.argv[0]+' <file> <variable> <LSM_file> <first record> <last record>'; sys.exit(0)
cf_in = sys.argv[1] ; cv_in=sys.argv[2] ; cf_lsm=sys.argv[3]

jt0 = int(sys.argv[4])
jtN = int(sys.argv[5])

roffset = 0.

if cv_in == 'curl_ssu':
    cfield = 'RV(SSU)'
    tmin=-40. ;  tmax=40.   ;  dtemp = 5.
    cpal_fld = 'ncview_blue_red'    
    cunit = r'$[10^{-6}s^{-1}]$'
    cb_jump = 1
    
if cv_in == 'mod_ssu':
    cfield = 'Surf. Current'
    tmin=0. ;  tmax=3.   ;  dtemp = 0.25
    cpal_fld = 'ncview_hotres'
    cunit = r'$[m/s]$'
    cb_jump = 1
    
if cv_in == 'sosstsst':
    cfield = 'SST'
    tmin=-2. ;  tmax=32.   ;  dtemp = 1.
    cpal_fld = 'ncview_nrl'    
    cunit = r'$^{\circ}C$'
    cb_jump = 2
    
elif cv_in == 'somxl010':
    cfield == 'MLD'
    tmin=50. ;  tmax=1500. ;  dtemp = 50.
    cpal_fld = 'viridis_r'
    

#bt.chck4f(cf_ice)

bt.chck4f(cf_in)
id_fld = Dataset(cf_in)
vtime = id_fld.variables['time_counter'][:]
id_fld.close()
Nt = len(vtime)

if jtN > Nt: print 'ERROR: your jtN > Nt !!!'; sys.exit(0)
if jtN == 0: jtN = Nt


bt.chck4f(cf_lsm)
id_lsm = Dataset(cf_lsm)
#LOLO for curl => F !!!
if lforce_mask: XMSK  = id_lsm.variables['fmask'][0,0,:,:]
XLON  =  id_lsm.variables['glamf'][0,:,:]
XLAT  =  id_lsm.variables['gphif'][0,:,:]
id_lsm.close()

(nj,ni) = nmp.shape(XLON)

#pmsk = nmp.ma.masked_where(XMSK[:,:] > 0.2, XMSK[:,:]*0.+40.)

#idx_oce = nmp.where(XMSK[:,:] > 0.5)


XFLD    = nmp.zeros((nj,ni))
if nsts>1:
    XFLDt   = nmp.zeros((nj,ni))
    XFLDtp1 = nmp.zeros((nj,ni))


params = { 'font.family':'Ubuntu',
           'font.size':       int(12),
           'legend.fontsize': int(12),
           'xtick.labelsize': int(12),
           'ytick.labelsize': int(12),
           'axes.labelsize':  int(12) }
mpl.rcParams.update(params)
cfont_clb   = { 'fontname':'Arial', 'fontweight':'normal', 'fontsize':13, 'color':'white' }
cfont_title = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':18, 'color':'white' }
cfont_mail  = { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':10, 'color':'0.5'}


# Colormaps for fields:
pal_fld = bcm.chose_colmap(cpal_fld)
norm_fld = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)

#pal_ice = bcm.chose_colmap(cpal_ice)
#norm_ice = colors.Normalize(vmin = rmin_ice, vmax = 1, clip = False)

pal_lsm = bcm.chose_colmap('blk')
norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)

vc_fld = nmp.arange(tmin, tmax + dtemp, dtemp)

print ''

rh = 7.5

jrot = -1 + jt0*nsts

print '\n jt0, jtN =>', jt0, jtN, '\n'

for jt in range(jt0,jtN):

    ct = '%3.3i'%(jt+1)

    cd = str(datetime.datetime.strptime('1990 '+ct, '%Y %j'))
    cdate = cd[:10] ; print ' *** Date :', cdate

    print "\n *** Reading record #"+str(ct)+" of "+cv_in+" in "+cf_in
    id_fld = Dataset(cf_in)

    if nsts==1:
        XFLD = id_fld.variables[cv_in][jt,:,:]
    else:
        if jt > jt0:
            XFLDt[:,:]  = XFLDtp1[:,:]
        else:
            XFLDt = id_fld.variables[cv_in][jt,:,:]
        XFLDtp1  = id_fld.variables[cv_in][jt+1,:,:]

    id_fld.close()
    print "  => done!"

    # sub time - step
    if nsts>1: dFdt = (XFLDtp1[:,:] - XFLDt[:,:])/float(nsts)
    
    for js in range(nsts):

        jrot = jrot+1
        
        # Linear interpolation:     # yN = y1 + (tN-t1)*(y2-y1)/(t2-t1)
        if nsts>1: XFLD[:,:] = XFLDt[:,:] + float(js)*dFdt
        
        rot = (long_start + (0.5/float(nsts)*float(jrot)))%360.
        rot = -rot

        cfig = 'figs/'+cv_in+'_NEMO'+'_d'+ct+'_'+str(js)+'.'+fig_type    
        
        print ' *** reference longitude =', rot

    	fig = plt.figure(num = 1, figsize=(rh,1.167*rh), dpi=None, facecolor='b', edgecolor='k')
    	ax  = plt.axes([0.005, 0.05, 0.99, 0.99], axisbg = '0.35')

        plt.title('Ocean (NEMO@ORCA12 + IFS@'+CTATM+'): '+cfield+', '+cdate, **cfont_title)

        print ' *** Creating new projection'
        carte = Basemap(projection='ortho', lat_0=latitude, lon_0=rot, resolution='h')
        x0,y0 = carte(XLON,XLAT)

        if lforce_mask: XFLD = nmp.ma.masked_where(XMSK[:,:] < 0.5, XFLD[:,:])
        
        print ' *** Ploting on map...'
        cf = carte.pcolor(x0, y0, XFLD, cmap=pal_fld, norm=norm_fld)
        if lforce_mask: cc = carte.contour(x0, y0, XMSK, [ 0. ], colors='k', linewidths=0.75)

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
        plt.savefig(cfig, dpi=160, orientation='portrait', facecolor='k')
        print '  => '+cfig+' created!\n'
        plt.close(1)

        del cf, carte, x0, y0

    if nsts>1: del dFdt
