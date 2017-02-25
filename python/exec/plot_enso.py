#!/usr/bin/env python

# L. Brodeau, April 2011

import sys
from string import replace
import os
from netCDF4 import Dataset

import barakuda_tool as bt
import barakuda_plot as bp

venv_needed = {'NN_SST','FIG_FORM'}
vdic = bt.check_env_var(sys.argv[0], venv_needed)

if len(sys.argv) != 2:
    print 'Usage: '+sys.argv[0]+' <monthly_SST_time_series.nc>'
    sys.exit(0)
cf_in  = sys.argv[1]

cname = replace(os.path.basename(cf_in), '.nc', '')

ff = vdic['FIG_FORM'] ; # format for figures (usually "png" or "svg")

bt.chck4f(cf_in)
id_in = Dataset(cf_in)
vtime = id_in.variables['time'][:] ; nbm = len(vtime)
vsst  = id_in.variables[vdic['NN_SST']][:]
id_in.close()

Nt = len(vsst)

ittic = bt.iaxe_tick(Nt/12)

bp.plot("enso")( vtime, vsst, cfignm=cname, dt=ittic, cfig_type=ff )
