# My own colormaps...

# Check http://matplotlib.org/examples/color/colormaps_reference.html !!!

# Colors:  http://www.pitt.edu/~nisg/cis/web/cgi/rgb.html

# Last Updated: L. Brodeau, January 2013

import sys
import numpy as nmp

# List of Barakuda home-made colormaps:
list_barakuda = [ 'blk', 'land', 'cb1', 'eke', 'bathy', 'mld', 'jetblanc', 'amoc',
                  'sst1', 'sst2', 'sst3', 'ice', 'blanc', 'rms',
                  'sigtr', 'bbr', 'bbr2', 'bbr0', 'bbr_cold', 'bbr_warm',
                  'cold0', 'warm0', 'graylb', 'graylb2', 'sigma', 'sigma0', 'mask' ]

# There is also NCVIEW colormaps, and the defauly Matplotlib colormaps...

l_debug = False

def chose_colmap( cname ):

    # 1st is it a ncview colormap ?
    if cname[:7] == 'ncview_':        
        M = ncview_cmap_to_array( cname )
        ColorMap = __build_colormap__(M)
        
        # Maybe a barakuda colormap ?
    elif cname in list_barakuda or ( cname[-2:] == '_r' and cname[:-2] in list_barakuda):
        if l_debug: print '\n *** Getting Barakuda colormap "'+cname+'" !'
        x = brkd_cmap(cname)
        ColorMap = x.clrmp()
        
    else:
        # Then it must be a Matplotlib colormap:
        from matplotlib.pylab import cm
        import matplotlib.pyplot as mp
        list = mp.colormaps()
        if cname in list:
            # Yes it is!
            if l_debug: print '\n *** Getting Matplotlib colormap "'+cname+'" !'
            fToCall = getattr(cm, cname)
            ColorMap = fToCall
        else:
            print 'ERROR: (chose_colmap of barakuda_colmap.py) do not know where to get colormap "'+cname+'" !'
            sys.exit(0)

    return ColorMap



def ncview_cmap_to_array( cname ):
    #
    #########################################################################################
    #
    # Get the NCVIEW colormap in the C header file and return an array ready for a colormap
    #         Author: L. Brodeau, 2017
    #
    # cname: "ncview_" +  name of the NCVIEW colormap as in the NCVIEW header colormap files  [string]
    #   example : 'ncview_rainbow'
    #      
    #
    # The environment variable 'DIR_NCVIEW_CMAP' must be set!
    #    => path to the directory containing the NCVIEW header colormap files  
    #    => ex: colormap 'rainbow' is defined in header file 'colormaps_rainbow.h''
    #
    ##########################################################################################
    
    import os
    import re

    dir_ncview_cmap = os.getenv('DIR_NCVIEW_CMAP')
    if dir_ncview_cmap is None:
        print(" ERROR => the {} environement variable is not set".format('DIR_NCVIEW_CMAP'))
        sys.exit(0)

    if cname[:7] != 'ncview_' : print ' ERROR: a ncview colormap should begin with "ncview_" !'; sys.exit(0)
    ncview_name = cname[7:]
    
    cf_ncview_cmap = dir_ncview_cmap+'/colormaps_'+ncview_name+'.h'
    if not os.path.exists(cf_ncview_cmap):
        print 'ERROR: NCVIEW colormap '+cf_ncview_cmap+' not found!' ; sys.exit(0)
    if l_debug: print '\n *** Getting NCVIEW colormap "'+ncview_name+'" from file "'+cf_ncview_cmap+'"'

    f = open(cf_ncview_cmap, 'r')
    cread_lines = f.readlines()
    f.close()

    lstarted = False
    vec = []
    for ll in cread_lines:

        ll = re.sub(r'\s', '', ll)
        ll = re.sub(r'{', ',', ll)
        ll = re.sub(r'}', ',', ll)
        ls = re.split(',',ll)
    
        if lstarted:
            for ve in ls[:-1]: vec.append(float(ve)) ; # [:-1] is to ommit the ',' at the end
    
        if ls[0] == 'staticintcmap_'+ncview_name+'[]=':
            lstarted = True        
            for ve in ls[1:-1]: vec.append(float(ve)) ; # [:-1] is to ommit the ',' at the end

    ctmp = []
    ii = 0
    while ii < len(vec):
        ctmp.append([vec[ii], vec[ii+1], vec[ii+2]])
        ii += 3

    MM = nmp.array(ctmp)/255.

    return MM








# ===== local ======


def __build_colormap__(MC, log_ctrl=0):

    import matplotlib.colors as mplc

    [ nc, n3 ] = nmp.shape(MC)

    # Make x vector :
    x =[]
    for i in range(nc): x.append(255.*float(i)/((nc-1)*255.0))
    x = nmp.array(x)
    if log_ctrl > 0: x = nmp.log(x + log_ctrl)
    rr = x[nc-1] ; x  = x/rr

    y =nmp.zeros(nc)
    for i in range(nc): y[i] = x[nc-1-i]

    x = 1 - y ; rr = x[nc-1] ; x  = x/rr

    vred  = [] ; vblue = [] ; vgreen = []

    for i in range(nc):
        vred.append  ([x[i],MC[i,0],MC[i,0]])
        vgreen.append([x[i],MC[i,1],MC[i,1]])
        vblue.append ([x[i],MC[i,2],MC[i,2]])

    cdict = {'red':vred, 'green':vgreen, 'blue':vblue}

    my_cm = mplc.LinearSegmentedColormap('my_colormap',cdict,256)

    return my_cm



#=======================================================================



class brkd_cmap:
    
    def __init__(self, name):
        self.name = name

    def clrmp(self):

        cname = self.name

        lrev = False        
        if cname[-2:] == '_r':
            lrev = True
            cname = cname[:-2]
        
        if cname == 'blk':
            M = nmp.array( [
                [ 0. , 0., 0. ], # black
                [ 0. , 0., 0. ]  # black
                ] )
            
        elif cname == 'land':
            M = nmp.array( [
                [ 0.4 , 0.4, 0.4 ],
                [ 0.4 , 0.4, 0.4 ]
                ] )

        elif cname == 'cb1':
            M = nmp.array( [
                [ 255,255,204  ],
                [ 161,218,180  ],
                [ 65,182,196  ],
                [ 44,127,184  ],
                [ 37,52,148  ]
                ] ) / 255.

        elif cname == 'eke':
            M = nmp.array( [
                [ 0.  , 0.0 , 0.2  ], # black
                [ 0.1 , 0.5 , 1.0  ], # blue
                [ 0.2 , 1.0 , 0.0  ], # green
                [ 1.  , 1.0 , 0.0  ], # yellow
                [ 1.  , 0.0 , 0.0  ], # red
                [0.2  , 0.27, 0.07 ] # brown
                ] )

        elif cname == 'bathy':
            M = nmp.array( [
                [ 0.0 , 0.0 , 0.4 ], # dark blue
                [ 0.1 , 0.5 , 1.0  ], # blue
                [ 0.2 , 1.0 , 0.0  ], # green
                [ 1.  , 1.0 , 0.0  ], # yellow
                [ 1.  , 0.0 , 0.0  ], # red
                [0.2  , 0.27, 0.07 ] # brown
                ] )

        elif cname == 'mld':
            M = nmp.array( [
                [ 1.0 , 1.0 , 1.0 ], # white
                [ 0.1 , 0.5 , 1.0 ], # light blue
                [ 0.13, 0.54, 0.13], # dark green
                [ 0.2 , 1.0 , 0.0 ], # light green
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 0.2 , 0.3 , 0.1 ] # dark redish brown
                ] )
        
        elif cname == 'jetblanc':
            M = nmp.array( [
                [ 0.6 , 0.0 , 0.8 ], # violet
                [ 0.0 , 0.0 , 0.4 ], # dark blue
                [ 0.1 , 0.5 , 1.0 ], # light blue
                [ 1.0 , 1.0 , 1.0 ], # white
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 0.2 , 0.3 , 0.1 ]  # dark redish brown
                ] )
        
        elif cname == 'amoc':
            M = nmp.array( [
                [ 0.4 , 0.0 , 0.6 ], # violet
                [ 1.0 , 1.0 , 1.0 ], # white
                [ 1.0 , 1.0 , 1.0 ], # white
                [0.68 , 0.98, 0.98], # light blue
                [ 0.0 , 0.0 , 0.95], # dark blue
                [ 0.2 , 1.0 , 0.0 ], # green
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 0.2 , 0.3 , 0.1 ] # dark read
                ] )
        
        elif cname == 'sst1':
            M = nmp.array( [
                [ 1.0 , 1.0 , 1.0 ], # white
                [ 0.4 , 0.0 , 0.6 ], # violet
                [ 0. , 0.2 , 0.99], # dark blue
                [0.68 , 0.98, 0.98], # light blue
                [ 0.13, 0.54, 0.13], # dark green
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 0.2 , 0.3 , 0.1 ] # dark read
                ] )
            
        elif cname == 'sst2':
            M = nmp.array( [
                [ 1.0 , 1.0 , 1.0 ], # white
                [ 0.4 , 0.0 , 0.6 ], # violet
                [ 0.0 , 0.0 , 0.95], # dark blue
                [0.68 , 0.98, 0.98], # light blue
                [46./255., 203./255., 35./255.], # green
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 0.2 , 0.3 , 0.1 ] # dark read
                ] )

        elif cname == 'sst3':
            M = nmp.array( [
                [ 0.4 , 0.0 , 0.6 ], # violet
                [ 0. , 0.2 , 0.99], # dark blue
                [0.68 , 0.98, 0.98], # light blue
                [ 0.13, 0.54, 0.13], # dark green
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 0.2 , 0.3 , 0.1 ] # dark read
                ] )
            
        elif cname == 'ice':
            M = nmp.array( [
                [  0. , 0.  , 0.3 ], # dark blue
                [ 0.6 , 0.6 , 0.8 ], # light grey
                [ 0.95 , 0.95 , 0.95 ],  # white
                [ 1.0 , 1.0 , 1.0 ]  # white
                ] )
        
        elif cname == 'blanc':
            M = nmp.array( [
                [ 1.0 , 1.0 , 1.0 ],  # white
                [ 1.0 , 1.0 , 1.0 ]  # white
                ] )
            
        elif cname == 'rms':
            M = nmp.array( [
                [ 1.0 , 1.0 , 1.0 ],
                [ 0.1 , 0.5 , 1.0 ],
                [ 0.2 , 1.0 , 0.0 ],
                [ 1.0 , 1.0 , 0.0 ],
                [ 1.0 , 0.0 , 0.0 ],
                [ 0.2 , 0.3 , 0.1 ]
                ] )

        elif cname == 'sigtr':
            M = nmp.array( [
                [ 1.0 , 1.0 , 1.0 ], # white
                [ 0.0 , 0.8 , 1.0 ], #light blue
                [ 0.1 , 0.5 , 1.0 ], #light blue
                [ 0.0 , 0.0 , 0.4 ], # blue
                [ 0.0 , 0.4 , 0.0 ], # dark green
                [ 0.1 , 1.0 , 0.0 ], # green
                [ 0.4 , 1.0 , 0.0 ], # vert pomme
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 1.0 , 0.4 , 0.0 ], # orange
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 0.6 , 0.0 , 0.0 ], # red
                [ 0.2 , 0.3 , 0.1 ]  # dark red
                ] )
            
        elif cname == 'bbr':
            M = nmp.array( [
                [ 0.  , 0. , 0.2 ],
                [ 0.  , 0. , 1.  ],
                [ 1.  , 1. , 1.  ],
                [ 1.  , 1. , 1.  ],
                [ 1.  , 0. , 0.  ],
                [ 0.6 , 0. , 0.  ]
                ] )
            
        elif cname == 'bbr2':
            M = nmp.array( [
                [ 0.  , 1. , 1.  ], # cyan
                [ 0.  , 0. , 1.  ],
                [ 1.  , 1. , 1.  ],
                [ 1.  , 1. , 1.  ],
                [ 1.  , 0. , 0.  ],
                [ 1.  , 1. , 0.  ]  # jaune
                ] )
            
        elif cname == 'bbr0':
            M = nmp.array( [
                [ 0.  , 1. , 1.  ], # cyan
                [ 0.  , 0. , 1.  ],
                [ 1.  , 1. , 1.  ],
                [ 1.  , 0. , 0.  ],
                [ 1.  , 1. , 0.  ]  # jaune
                ] )
            
        elif cname == 'bbr_cold':
            M = nmp.array( [
                [ 0.  , 1. , 1.  ], # cyan
                [ 0.  , 0. , 1.  ],
                [ 19./255.  , 7./255. , 129./255  ], # dark blue
                [ .1  , .1 , .9  ],
                [ 1.  , 1. , 1.  ],
                [ 1.  , 1. , 1.  ],
                [ 0.7  , 0. , 0.  ] # Dark red
                ] )
            
        elif cname == 'bbr_warm':
            M = nmp.array( [
                [ 19./255.  , 7./255. , 129./255  ], # dark blue
                [ 1.  , 1. , 1.  ],
                [ 1.  , 1. , 1.  ],
                [ 0.9  , 0.1 , 0.1  ],
                [ 0.7  , 0. , 0.  ], # Dark red
                [ 1.  , 1. , 0.  ],  # jaune
                ] )
            
        elif cname == 'cold0':
            M = nmp.array( [
                [ 177./255.  , 250./255. , 122./255. ],   # greenish
                [ 0.  , 1. , 1.  ], # cyan
                [ 7./255.  , 11./255. , 122./255. ], # dark blue
                [ 0.  , 0. , 1.  ], # true blue
                [ 177./255.  , 189./255. , 250./255. ], # light blue
                [ 1.  , 1. , 1.  ],
                ] )
            
        elif cname == 'warm0':
            M = nmp.array( [
                [ 1.  , 1. , 1.  ],
                [ 255./255.  , 254./255. , 198./255.  ], # very light yellow
                [ 1.  , 1. , 0.  ],  # yellow
                [ 244./255.  , 78./255. , 255./255.  ], # pink
                [ 1.  , 0. , 0.  ], # true red
                [ 139./255.  , 5./255. , 5./255.  ] # dark red
                ] )
            
        elif cname == 'graylb':
            M = nmp.array( [
                [ 1.  , 1. , 1. ],
                [ 0.1  , 0.1 , 0.1 ]
                ] )
            
        elif cname == 'graylb2':
            M = nmp.array( [
                [ 0.6  , 0.6 , 0.6 ],
                [ 1.  , 1. , 1. ]
                ] )
            
        elif cname == 'sigma':
            M = nmp.array( [
                [ 1.0 , 1.0 , 1.0 ], # white
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 0.2 , 1.0 , 0.0 ], # green
                [ 0.1 , 0.5 , 1.0 ], # light blue
                [ 0.0 , 0.0 , 0.4 ], # dark blue
                [ 0.6 , 0.0 , 0.8 ]  # violet
                ] )

        elif cname == 'sigma0':
            M = nmp.array( [
                [ 0.2 , 0.3 , 0.1 ], # dark redish brown
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 0.2 , 1.0 , 0.0 ], # green
                [ 0.1 , 0.5 , 1.0 ], # light blue
                [ 0.0 , 0.0 , 0.4 ], # dark blue
                [ 0.6 , 0.0 , 0.8 ], # violet
                [ 1.0 , 1.0 , 1.0 ]  # white
                ] )
            
        elif cname == 'mask':
            M = nmp.array( [
                [ 0.5 , 0.5 , 0.5 ], # gray
                [ 0.5 , 0.5 , 0.5 ]  # gray
                ] )

        else:
            print 'ERROR: (''barakuda_colmap.py) => unknown "barakuda" colormap: '+cname
            sys.exit(0)

        if lrev:
            # reverse colormap:
            my_cmap = __build_colormap__(M[::-1,:])
        else:
            my_cmap = __build_colormap__(M)
             
        return my_cmap

#=======================================================================
