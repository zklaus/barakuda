#!/usr/bin/env python

#       B a r a K u d a
#
#
#
#       L. Brodeau, 2017
#

import sys
import os
#import numpy as nmp
import barakuda_tool as bt
#import barakuda_ncio as bn
#import barakuda_orca as bo
#import barakuda_plot as bp


csn = sys.argv[0]

venv_needed = {'ORCA','EXP'} ; #,'TRANSPORT_SECTION_FILE'}

vdic = bt.check_env_var(csn, venv_needed)

CONFEXP = vdic['ORCA']+'-'+vdic['EXP']






















narg = len(sys.argv)
if narg != 2:
    print 'Usage: {} <ASCII section file>'.format(csn)
    sys.exit(0)
cf_in = sys.argv[1]


bt.chck4f(cf_in, script_name=csn)                


list_sections = bt.get_sections_from_file(cf_in)



print '\n\n Section:\n', list_sections[:]
























print ''+csn+' done...\n'






