#!/usr/bin/env python

# L. Brodeau, February 2017

# Plot the Atlantic Multidecadal Oscillation out of a time-series of SST
# averaged over the North Atlantic (0-70N)

import sys
from string import replace
import os
from netCDF4 import Dataset
import numpy as nmp

import barakuda_tool as bt
import barakuda_plot as bp
import barakuda_stat as bs

n_run_mean = 11
#n_run_mean = 21

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
vt_m   = id_in.variables['time'][:]
vsst_m = id_in.variables[cv_in][:]
id_in.close()

Nm = len(vsst_m)

# Annual averaging first:
vt, vsst = bt.monthly_2_annual(vt_m, vsst_m)

Nt = len(vt)

n1  = (n_run_mean-1)/2
n2  = -n1 - 1

vtmp = nmp.zeros(Nt)
vtime = nmp.zeros(Nt-n_run_mean)
xplot = nmp.zeros((Nt-n_run_mean,4)) ; # Nt-n_run_mean because X month-running mean

vtime[:] = vt[n1:n2]
if n_run_mean == 11: vtmp = bs.running_mean_11(vsst) ; # 11-month running mean
if n_run_mean == 21: vtmp = bs.running_mean_21(vsst) ; # 21-month running mean
xplot[:,0] = vsst[n1:n2]
xplot[:,1] = vtmp[n1:n2]

(za,zb) = bs.least_sqr_line(vtime[:], xplot[:,1]) ; # Least-square linear trend
xplot[:,2] = za*vtime[:] + zb

xplot[:,3] = xplot[:,1] - xplot[:,2] ; # anomaly for 11-month running mean

ittic = bt.iaxe_tick(Nm/12)

bp.plot("oscillation_index")( vtime, xplot[:,3], ymax=0.3, dy=0.05,
                              tmin=vt_m[0], tmax=vt_m[-1], dt=ittic,
                              cfignm=cname, cfig_type=fig_form,
                              cyunit=r'SST anomaly ($^{\circ}$C)',
                              ctitle='Atlantic Multidecadal Oscillation ('+str(n_run_mean)+'-year running mean)' )

