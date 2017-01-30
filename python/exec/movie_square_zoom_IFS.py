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
import gc

import barakuda_colmap as bcm

import barakuda_tool as bt

gc.collect()

# South Greenland:
#i1 = 412; i2 =486
#j1 = 22 ; j2 = 56

# NAtl:
#i1 = 385 ; i2= 540
#j1 =   6 ; j2 = 84

#Global T255:
#i1 = 0 ; i2 =511
#j1 = 0 ; j2 = 255

#Global T1279:
i1 = 0 ; i2 = 2560
j1 = 0 ; j2 = 1280

year_ref_ini = 1990

fig_type='png'

narg = len(sys.argv)
if narg < 4: print 'Usage: '+sys.argv[0]+' <file> <variable> <LSM_file>'; sys.exit(0)
cf_in = sys.argv[1] ; cv_in=sys.argv[2] ; cf_lsm=sys.argv[3]

lsst = False ; lshf = False 
if cv_in == 'SSTK': lsst = True
if cv_in == 'SNHF': lshf = True


if lsst:
    tmin=-20.  ;  tmax=12. ;  dt = 1.
    cpal = 'sstnw'
    cfield = 'SST'
    
if lshf:
    tmin=-1200. ;  tmax=400. ;  dt = 25.
    #cpal = 'rainbow'
    cpal = 'nrl'
    cfield = 'Net Heat Flux'


clsm  = 'LSM'


imax=i2+1

Ni = i2-i1
Nj = j2-j1


LSM = nmp.zeros((Nj,Ni), dtype=nmp.float)
XIN  = nmp.zeros((Nj,Ni))



bt.chck4f(cf_lsm)
id_lsm = Dataset(cf_lsm)

 #    i1      imax    i2

if i2 >= imax:
    print ' i2 > imax !!! => ', i2, '>', imax
    Xall = id_lsm.variables[clsm][0,j1:j2,:]
    LSM[:,0:imax-i1] = Xall[:,i1:imax]
    ii=imax-i1
    LSM[:,ii:Ni] = Xall[:,0:i2-imax]
    del Xall
else:
    LSM[:,:]  = id_lsm.variables[clsm][0,j1:j2,i1:i2]
id_lsm.close()


[ nj , ni ] = nmp.shape(LSM)

idx_ocean = nmp.where(LSM[:,:] < 0.5)
LSM[idx_ocean] = nmp.nan
LSM = nmp.flipud(LSM)


params = { 'font.family':'Ubuntu',
           'font.size':       int(15),
           'legend.fontsize': int(15),
           'xtick.labelsize': int(15),
           'ytick.labelsize': int(15),
           'axes.labelsize':  int(15) }
mpl.rcParams.update(params)

# Pal_Sst:
pal_sst = bcm.ncview_colmap( cpal, '/home/Earth/lbrodeau/DEV/barakuda/src/ncview_colormaps' )
#pal_sst = bcm.chose_palette(cpal)
norm_sst = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)

pal_lsm  = bcm.chose_palette('blk')
norm_lsm = colors.Normalize(vmin = 0, vmax = 1, clip = False)


vc_sst = nmp.arange(tmin, tmax + dt, dt)

pfin = nmp.zeros((nj,ni))
    


bt.chck4f(cf_in)
id_in = Dataset(cf_in)
vtime = id_in.variables['time'][:]
id_in.close()
del id_in
Nt    = len(vtime)


for jt in range(Nt):

    print '\n *** Reading record # '+str(jt+1)+' of '+cv_in+' in '+cf_in
    id_in = Dataset(cf_in)
    if i2 >= imax:
        print ' i2 = ', i2
        Xall = id_in.variables[cv_in][jt,j1:j2,:]
        XIN[:,0:imax-i1] = Xall[:,i1:imax]
        ii=imax-i1
        XIN[:,ii:Ni] = Xall[:,0:i2-imax]
        del Xall
    else:
        XIN[:,:] = id_in.variables[cv_in][jt,j1:j2,i1:i2]
    id_in.close()
    del id_in

    if lsst: XIN[:,:] = XIN[:,:] - 273.15

    ct = '%3.3i'%(jt+1)

    cd = str(datetime.datetime.strptime(str(year_ref_ini)+' '+ct, '%Y %j'))
    cdate = cd[:10] ; print ' *** cdate :', cdate

    cfig = 'figs/'+cv_in+'_IFS'+'_d'+ct+'.'+fig_type

    fig = plt.figure(num = 1, figsize=(8.*float(Ni)/float(Nj)*0.8 , 8.), dpi=None, facecolor='w', edgecolor='k')
    #ax  = plt.axes([0.04, -0.06, 0.93, 1.02], axisbg = 'k')
    ax  = plt.axes([0.04, -0.04, 0.93, 1.02], axisbg = 'k')

    cf = plt.imshow(nmp.flipud(XIN), cmap = pal_sst, norm = norm_sst)

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

    del cf

    # Mask
    print ' LSM stuff...'    
    cm = plt.imshow(LSM, cmap = pal_lsm, norm = norm_lsm)
    print ' ... done!\n'

    cfont = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':22 }
    plt.title('IFS: '+cfield+', coupled ORCA12-T255, '+cdate, **cfont)

    
    plt.savefig(cfig, dpi=120, orientation='portrait', transparent=False)
    print cfig+' created!\n'
    plt.close(1)
                
    del fig, ax, clb, cm

    gc.collect()
    


