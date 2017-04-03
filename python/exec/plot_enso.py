#!/usr/bin/env python

# L. Brodeau, April 2011

import sys
from string import replace
import os
import numpy as nmp
from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_plot as bp
import barakuda_stat as bs


fig_form = os.getenv('FIG_FORM')
if fig_form is None: fig_form = 'png'

if len(sys.argv) != 3:
    print 'Usage: '+sys.argv[0]+' <monthly_SST_time_series.nc> <name_sst>'
    sys.exit(0)
cf_in  = sys.argv[1]
cv_in  = sys.argv[2]

bt.chck4f(cf_in)

cname = replace(os.path.basename(cf_in), '.nc', '')

id_in = Dataset(cf_in)
vt    = id_in.variables['time'][:]
vsst  = id_in.variables[cv_in][:]
id_in.close()

Nt = len(vsst)


vtmp = nmp.zeros(Nt)
vtime = nmp.zeros(Nt-5)
xplot = nmp.zeros((Nt-5,4)) ; # Nt-5 because 5-month-running mean

vtime[:] = vt[2:-3]
vtmp = bs.running_mean_5(vsst) ; # 5-month running mean
xplot[:,0] = vsst[2:-3]
xplot[:,1] = vtmp[2:-3]

(za,zb) = bs.least_sqr_line(vtime[:], xplot[:,1]) ; # Least-square linear trend
xplot[:,2] = za*vtime[:] + zb

xplot[:,3] = xplot[:,1] - xplot[:,2] ; # anomaly for 5-month running mean

ittic = bt.iaxe_tick(Nt/12)

bp.plot("oscillation_index")( vtime, xplot[:,3], ymax=2.1, dy=0.5, yplusminus=0.4, dt=ittic,
                              cfignm=cname, cfig_type=fig_form,
                              cyunit=r'SST anomaly ($^{\circ}$C)',
                              ctitle='ENSO (over Nino box 3.4)' )
