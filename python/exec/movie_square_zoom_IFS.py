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

year_ref_ini = 1990

#CTATM = 'T255'
CTATM = 'T1279'

# NCVIEW colormaps?
dir_ncview_cmap = os.getenv('DIR_NCVIEW_CMAP')
if dir_ncview_cmap is None:
    print(" ERROR => the {} environement variable is not set".format('DIR_NCVIEW_CMAP'))
    sys.exit(0)

if CTATM == 'T255':
    # South Greenland:
    #i1 = 412; i2 =486
    #j1 = 22 ; j2 = 56    
    # NAtl:
    i1 = 385 ; i2= 540
    j1 =   6 ; j2 = 84
    #Global T255:
    #i1 = 0 ; i2 =511
    #j1 = 0 ; j2 = 255

elif CTATM == 'T1279':
    #
    #Global:
    #i1 = 0 ; i2 = 2559+1
    #j1 = 0 ; j2 = 1279+1
    # Natl:
    #i1 = 1849 ; i2 = 2525
    #j1 = 97   ; j2 = 508
    i1 = 1960 ; i2 = 2550; #2680
    j1 = 0   ; j2 = 519

else:
    print 'UNKNOW ATMOSPHERE RESOLUTION!'; sys.exit(0)


fig_type='png'

narg = len(sys.argv)
if narg < 4: print 'Usage: '+sys.argv[0]+' <file> <variable> <LSM_file>'; sys.exit(0)
cf_in = sys.argv[1] ; cv_in=sys.argv[2] ; cf_lsm=sys.argv[3]

lsst = False ; lshf = False
if cv_in == 'T2M':  lt2m = True
if cv_in == 'SSTK': lsst = True
if cv_in == 'SNHF': lshf = True


if lt2m:
    tmin=-16.  ;  tmax=28. ;  dt = 1.
    #cpal = 'nrl'
    #cpal = 'jaisnd'
    #cpal = '3gauss'
    #cpal = 'rainbow2_cmyk'
    cpal = 'rainbow'
    #cpal = 'rnb2'
    #cpal = 'jaisnc'
    #cpal = 'jaisnb'
    cfield = 'T2M'
    cunit = r'$^{\circ}C$'
    cb_jump = 2
    
if lsst:
    tmin=-20.  ;  tmax=12. ;  dt = 1.
    cpal = 'sstnw'
    cfield = 'SST'
    cunit = r'$Boo$'
    cb_jump = 2
    
if lshf:
    tmin=-1200. ;  tmax=400. ;  dt = 25.
    #cpal = 'rainbow'
    cpal = 'nrl'
    cfield = 'Net Heat Flux'
    cunit = r'$W/m^2$'
    cb_jump = 4


clsm  = 'LSM'


# Need to know dimension:
bt.chck4f(cf_lsm)
id_lsm = Dataset(cf_lsm)
vlon = id_lsm.variables['lon'][:]
vlat = id_lsm.variables['lat'][:]
id_lsm.close()

Ni0 = len(vlon)
Nj0 = len(vlat)

print '\n Dimension of global domain:', Ni0, Nj0


imax=Ni0+1
Ni = i2-i1
Nj = j2-j1


LSM = nmp.zeros((Nj,Ni), dtype=nmp.float)
XIN  = nmp.zeros((Nj,Ni))


id_lsm = Dataset(cf_lsm)
if i2 >= imax:
    print ' i2 > imax !!! => ', i2, '>', imax
    Xall = id_lsm.variables[clsm][0,j1:j2,:]
    LSM[:,0:imax-i1] = Xall[:,i1-1:imax]
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
           'font.size':       int(12),
           'legend.fontsize': int(12),
           'xtick.labelsize': int(12),
           'ytick.labelsize': int(12),
           'axes.labelsize':  int(12) }
mpl.rcParams.update(params)
cfont_clb   = { 'fontname':'Arial', 'fontweight':'normal', 'fontsize':13 }
cfont_title = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':18 }
cfont_mail  = { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':9, 'color':'0.5' }



# Pal_Sst:
pal_fld = bcm.ncview_colmap( cpal, dir_ncview_cmap )
#pal_fld = bcm.chose_palette(cpal)
norm_fld = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)

pal_lsm  = bcm.chose_palette('blk')
norm_lsm = colors.Normalize(vmin = 0, vmax = 1, clip = False)


vc_fld = nmp.arange(tmin, tmax + dt, dt)

pfin = nmp.zeros((nj,ni))
    


bt.chck4f(cf_in)
id_in = Dataset(cf_in)
vtime = id_in.variables['time'][:]
id_in.close()
del id_in
Nt    = len(vtime)


# Size of the figure:
rat_Nj_Ni = float(Nj)/float(Ni) + 0.12
rh =  7.5
rw = rh/rat_Nj_Ni
FSZ = ( rw  , rh )
rcorr = rat_Nj_Ni/(float(Nj0)/float(Ni0))
print '  rcorr => ', rcorr


for jt in range(Nt):

    print '\n *** Reading record # '+str(jt+1)+' of '+cv_in+' in '+cf_in
    id_in = Dataset(cf_in)
    if i2 >= imax:
        print ' i2 = ', i2
        Xall = id_in.variables[cv_in][jt,j1:j2,:]
        XIN[:,0:imax-i1] = Xall[:,i1-1:imax]
        ii=imax-i1
        XIN[:,ii:Ni] = Xall[:,0:i2-imax]
        del Xall
    else:
        XIN[:,:] = id_in.variables[cv_in][jt,j1:j2,i1:i2]
    id_in.close()
    del id_in

    if lsst or lt2m: XIN[:,:] = XIN[:,:] - 273.15

    ct = '%3.3i'%(jt+1)

    cd = str(datetime.datetime.strptime(str(year_ref_ini)+' '+ct, '%Y %j'))
    cdate = cd[:10] ; print ' *** cdate :', cdate

    cfig = 'figs/'+cv_in+'_IFS'+'_d'+ct+'.'+fig_type

    
    fig = plt.figure(num = 1, figsize=FSZ, dpi=None, facecolor='w', edgecolor='k')

    ax  = plt.axes([0.055, 0.05, 0.9, 1.], axisbg = 'k')

    cf = plt.imshow(nmp.flipud(XIN), cmap = pal_fld, norm = norm_fld)

    plt.axis([ 0, ni, 0, nj])

    # Mask
    print ' LSM stuff...'    
    cm = plt.imshow(LSM, cmap = pal_lsm, norm = norm_lsm)

    plt.title('IFS: '+cfield+', coupled ORCA12-'+CTATM+', '+cdate, **cfont_title)







    ax2 = plt.axes([0.04, 0.08, 0.93, 0.025])
    clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend='both')
    #clb = plt.colorbar(cf, ticks=vc_fld, orientation='horizontal', drawedges=False, pad=0.07, shrink=1., aspect=40)
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


    ax.annotate('laurent.brodeau@bsc.es', xy=(1, 4), xytext=(480, -85), **cfont_mail)
    
    plt.savefig(cfig, dpi=160, orientation='portrait', transparent=False)
    print cfig+' created!\n'
    plt.close(1)
                
    del fig, ax, clb, cm

    


