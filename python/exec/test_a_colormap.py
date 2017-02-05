#!/usr/bin/env python

#       B a r a K u d a
#
#    L. Brodeau, 2017
#

import sys
import numpy as nmp

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mplc

import barakuda_colmap as bcm



narg = len(sys.argv)
if narg < 3:
    print 'Usage: '+sys.argv[0]+' <Colormap> <Min val> <Max val> <increment>'
    print '      ex: $ test_a_colormap.py rainbow    -10 50 5'
    print '          $ test_a_colormap.py ncview_nrl  0  20 1'
    sys.exit(0)

cpal_fld      = sys.argv[1]
tmin    = float(sys.argv[2])
tmax    = float(sys.argv[3])
rdt     = float(sys.argv[4])

vc_fld = nmp.arange(tmin, tmax + rdt, rdt)


#for cpal_fld in bcm.list_cmap_barakuda :

for cpal_fld in [ cpal_fld ] :

    cpal_fld = cpal_fld

    cfig = 'show_colmap_'+cpal_fld+'.png'

    pal_fld = bcm.chose_colmap(cpal_fld)

    nrm_fld = mplc.Normalize(vmin = tmin, vmax = tmax, clip = False)

    fig = plt.figure(num = 1, figsize=(10,2), dpi=None, facecolor='w', edgecolor='k')

    ax  = plt.axes([0.01, 0.12, 0.98, 0.8])

    clb = mpl.colorbar.ColorbarBase(ax, ticks=vc_fld, cmap=pal_fld, norm=nrm_fld,
                                    orientation='horizontal', extend='both')

    plt.savefig(cfig, dpi=120, orientation='portrait', transparent=False)
    plt.close(1)

    print '\n Check how "'+cpal_fld+'" looks like in: '+cfig+' !\n'


#cb_labs = [] ; cpt = 0
#for rr in vc_fld:
#    if cpt % cb_jump == 0:
#        cb_labs.append(str(int(rr)))
#    else:
#        cb_labs.append(' ')
#    cpt = cpt + 1
#clb.ax.set_xticklabels(cb_labs)
#clb.set_label(cunit, **cfont_clb)
