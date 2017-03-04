#!/usr/bin/env python
#
##############################################################
#       B a r a K u d a
#
#        Generate netcdf files of cross-sections
#
#       L. Brodeau, 2016
##############################################################

import sys
import numpy as nmp
from netCDF4 import Dataset

import barakuda_orca as bo
import barakuda_tool as bt
import barakuda_ncio as bnc

venv_needed = {'ORCA','EXP','DIAG_D','i_do_sect','TS_SECTION_FILE','MM_FILE','NN_T','NN_S'}
vdic = bt.check_env_var(sys.argv[0], venv_needed)


i_do_sect = int(vdic['i_do_sect'])

if i_do_sect != 1: print 'ERROR: sys.argv[0] => why are we here when i_do_sect != 1 ???'; sys.exit(0)

f_sections = vdic['TS_SECTION_FILE']

CONFEXP = vdic['ORCA']+'-'+vdic['EXP']

cnexec = sys.argv[0]

na = len(sys.argv)
if na != 3:
    print 'Usage : '+cnexec+' <EXP_grid_T.nc> <year>'
    sys.exit(0)

cf_in  = sys.argv[1]
cyear  = sys.argv[2] ; jyear = int(cyear); cyear = '%4.4i'%jyear


cv_t = vdic['NN_T']
cv_s = vdic['NN_S']

print 'Current year is '+cyear+' !\n'


bt.chck4f(vdic['MM_FILE'])
id_mm = Dataset(vdic['MM_FILE'])
rmsk = id_mm.variables['tmask'][0,:,:,:]
xlon = id_mm.variables['glamt'][0,:,:]
xlat = id_mm.variables['gphit'][0,:,:]
id_mm.close()

[ nk, nj, ni ] = rmsk.shape

bt.chck4f(cf_in)
id_in  = Dataset(cf_in)
vdepth = id_in.variables['deptht'][:]
XT     = id_in.variables[cv_t][:,:,:,:]
XS     = id_in.variables[cv_s][:,:,:,:]
id_in.close()



[ Nt, nk0, nj0, ni0 ] = XT.shape

if [ nk0, nj0, ni0 ] != [ nk, nj, ni ]: print 'ERROR: ssx_boxes.py => mask and field disagree in shape!'; sys.exit(0)

print 'Nt, nk, nj, ni =', Nt, nk, nj, ni


# Masking:
for jt in range(Nt):
    XT[jt,:,:,:] = rmsk[:,:,:]*XT[jt,:,:,:] + (1. - rmsk[:,:,:])*-9999.
    XS[jt,:,:,:] = rmsk[:,:,:]*XS[jt,:,:,:] + (1. - rmsk[:,:,:])*-9999.



vtime = nmp.zeros(Nt)
for jt in range(Nt): vtime[jt] = float(jyear) + (float(jt) + 0.5)/float(Nt)


# Getting sections:
vboxes, vlon1, vlat1, vlon2, vlat2 = bt.read_coor(f_sections, ctype='float', lTS_bounds=False)


js = -1
for csname in vboxes:

    js = js + 1
    
    print'\n *** '+sys.argv[0]+': treating section '+csname
    
    ( i1, i2, j1, j2 ) = bo.transect_zon_or_med(vlon1[js], vlon2[js], vlat1[js], vlat2[js], xlon, xlat)

    print csname+' :'
    print   '(lon1, lon2, lat1, lat2) =', vlon1[js], vlon2[js], vlat1[js], vlat2[js]
    print   ' => i1, i2, j1, j2 =', i1, i2, j1, j2
    print ''

    if i1 > i2: print 'ERROR: cross_sections.py => i1 > i2 !'; sys.exit(0)
    if j1 > j2: print 'ERROR: cross_sections.py => j1 > j2 !'; sys.exit(0)


    if i1 == i2:
        print 'Meridional section!'
        caxis = 'y' ; cxn = 'lat'         
        vaxis = xlat[j1:j2,i1]
        imsk  = rmsk[:,j1:j2,i1]
        ZT = XT[:,:,j1:j2,i1]
        ZS = XS[:,:,j1:j2,i1]

    if j1 == j2:
        print 'Zonal section!'
        caxis = 'x'; cxn = 'lon' 
        vx = xlon[j1,i1:i2] ; vaxis = nmp.zeros(len(vx)) ; vaxis[:] = vx[:]
        ivf = nmp.where(vx>180); vaxis[ivf] = vx[ivf] - 360.
        imsk  = rmsk[:,j1,i1:i2]
        ZT = XT[:,:,j1,i1:i2]
        ZS = XS[:,:,j1,i1:i2]


    cf_out = vdic['DIAG_D']+'/TS_section_'+csname+'.nc'

    bnc.wrt_appnd_2dt_series(vaxis, -vdepth, vtime, ZT, cf_out, cv_t,
                             missing_val=-9999.,
                             cxdnm=cxn, cydnm='depth', cxvnm=cxn, cyvnm='depth',
                             cu_t='year', cu_d='deg.C', cln_d='Potential temperature',
                             xd2=ZS, cvar2=cv_s, cln_d2='Salinity', cun2='PSU')
