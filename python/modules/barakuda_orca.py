
import sys
import numpy as nmp



# Definition of the boxes used to average MLD on
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  => used into mean.py and plot_time_series.py and compare_time_series.py...
cdgS=r'$^{\circ}$S'
cdgN=r'$^{\circ}$N'

cname_mld_boxes = [ 'Glob'  ,    'Sof40S'        ,      'Sof55S'      ,      'Nof40N'      , 'NAtl'                    , '40S40N'                   ]
clgnm_mld_boxes = [ 'Global', 'South of 40'+cdgS , 'South of 55'+cdgS , 'North of 40'+cdgN , 'Northern North Atlantic' , '40$^{\circ}$S - 40'+cdgN  ]
                                                                                                           
r_lon_p1_mld    = [   -999. ,     -999.          ,        -999.       ,        -999.       ,  -75.                     ,     -999.                  ]
r_lon_p2_mld    = [   -999. ,     -999.          ,        -999.       ,        -999.       ,   15.                     ,     -999.                  ]
                                                                                                
r_lat_p1_mld    = [   -999. ,     -999.          ,        -999        ,         40.        ,   50.                     ,      -40.                  ]
r_lat_p2_mld    = [   -999. ,      -40.          ,        -55.        ,        -999.       ,   68.                     ,       40.                  ]






def get_basin_info( cf_bm ):    
    from netCDF4 import Dataset
    l_b_names = [] ; l_b_lgnms = []
    l_b_names.append(u'GLO') ; l_b_lgnms.append(u'Global Ocean')
    id_bm = Dataset(cf_bm)
    list_var = id_bm.variables.keys()
    for cv in list_var:
        if cv[:5] == 'tmask':
            l_b_names.append(cv[5:])
            l_b_lgnms.append(id_bm.variables[cv].long_name)
    id_bm.close()
    return l_b_names, l_b_lgnms


def lon_reorg_orca(ZZ, corca, ilon_ext):
    #
    #
    # IN:
    # ===
    # ZZ       : array to  longitude, 3D, 2D field or 1D longitude vector
    # corca    : ORCA conf
    # ilon_ext : longitude extention in degrees
    #
    # OUT:
    # ====
    # ZZx     : re-organized array, the x dimension is now nx-2
    #
    import barakuda_tool as bt
    #
    # jx_junc : ji when lon become positive!
    if   corca[:5] == 'ORCA2':
        jx_junc = 141
    elif corca[:5] == 'ORCA1' and corca[5] != '2':
        jx_junc = 288
    elif corca[:6] == 'eORCA1' and corca[6] != '2':
        jx_junc = 288
    elif corca[:7] == 'ORCA025':
        jx_junc = 1150
    else:
        print 'ERROR: lon_reorg_orca.barakuda_orca => '+corca+' not supported yet!'; sys.exit(0)

    jx_oo = 2  # orca longitude overlap...

    vdim = ZZ.shape

    ndim = len(vdim)

    if ndim < 1 or ndim > 4:
        print 'util_orca.lon_reorg_orca: ERROR we only treat 1D, 2D, 3D or 4D arrays...'

    if ndim == 4:
        #
        [ nr, nz , ny , nx ] = vdim ;     nx0 = nx - jx_oo
        ZZx  = nmp.zeros((nr, nz, ny, nx0))
        ZZx_ext  = nmp.zeros((nr, nz, ny, (nx0+ilon_ext)))
        #
        for jx in range(jx_junc,nx):
            ZZx[:,:,:,jx-jx_junc] = ZZ[:,:,:,jx]
        for jx in range(jx_oo,jx_junc):
            ZZx[:,:,:,jx+(nx-jx_junc)-jx_oo] = ZZ[:,:,:,jx]
        #
        if ilon_ext == 0: ZZx_ext[:,:,:,:] = ZZx[:,:,:,:]
    #
    #
    if ndim == 3:
        #
        [ nz , ny , nx ] = vdim ;     nx0 = nx - jx_oo
        #print 'nx, ny, nz = ', nx, ny, nz
        #
        ZZx  = nmp.zeros(nx0*ny*nz) ;  ZZx.shape = [nz, ny, nx0]
        ZZx_ext  = nmp.zeros((nx0+ilon_ext)*ny*nz) ;  ZZx_ext.shape = [nz, ny, (nx0+ilon_ext)]
        #
        for jx in range(jx_junc,nx):
            ZZx[:,:,jx-jx_junc] = ZZ[:,:,jx]
        for jx in range(jx_oo,jx_junc):
            ZZx[:,:,jx+(nx-jx_junc)-jx_oo] = ZZ[:,:,jx]
        #
        if ilon_ext == 0: ZZx_ext[:,:,:] = ZZx[:,:,:]
    #
    #
    if ndim == 2:
        #
        [ ny , nx ] = vdim ;     nx0 = nx - jx_oo
        #print 'nx, ny = ', nx, ny
        #
        ZZx  = nmp.zeros(nx0*ny) ;  ZZx.shape = [ny, nx0]
        ZZx_ext  = nmp.zeros((nx0+ilon_ext)*ny) ;  ZZx_ext.shape = [ny, (nx0+ilon_ext)]
        #
        for jx in range(jx_junc,nx):
            ZZx[:,jx-jx_junc] = ZZ[:,jx]
        for jx in range(jx_oo,jx_junc):
            ZZx[:,jx+(nx-jx_junc)-jx_oo] = ZZ[:,jx]
        #
        if ilon_ext == 0: ZZx_ext[:,:] = ZZx[:,:]
    #
    #
    if ndim == 1:
        #
        [ nx ] = vdim ;     nx0 = nx - jx_oo
        #print 'nx = ', nx
        #
        ZZx  = nmp.zeros(nx0) ;  ZZx.shape = [nx0]
        ZZx_ext  = nmp.zeros(nx0+ilon_ext) ;  ZZx_ext.shape = [nx0+ilon_ext]
        #
        for jx in range(jx_junc,nx):
            ZZx[jx-jx_junc] = ZZ[jx]
            #print jx-jx_junc, 'prend', jx, '    ', vlon[jx]
            #
        #print ''
        for jx in range(jx_oo,jx_junc):
            ZZx[jx+(nx-jx_junc)-jx_oo] = ZZ[jx]
            #print jx+(nx-jx_junc)-jx_oo, 'prend', jx, '    ', vlon[jx]
        #
        if ilon_ext == 0: ZZx_ext[:] = ZZx[:]
        #iwa = nmp.where(vlon0 < 0.) ; vlon0[iwa] = vlon0[iwa] + 360.
        #
        #
        #
        # Now longitude extenstion:
    if ilon_ext > 0: ZZx_ext = bt.extend_domain(ZZx, ilon_ext)
    del ZZx

    return ZZx_ext




def conf_run(ccr):
    #
    # Find the CONF from CONF-RUN and exit if CONF does not exist!
    #
    i = 0 ; conf = ''
    while i < len(ccr) and ccr[i] != '-' : conf = conf+ccr[i]; i=i+1
    #print 'conf =', conf, '\n'
    return conf


def info_run(ccr):
    #
    i = 0 ; j = 0 ; conf = '' ; case = '' ; cfrq = '' ; cyyy = ''
    #
    while i < len(ccr) and ccr[i] != '-' : conf = conf+ccr[i]; i=i+1
    i=i+1
    while i < len(ccr) and ccr[i] != '_' : case = case+ccr[i]; i=i+1
    i=i+1
    while i < len(ccr) and ccr[i] != '_' : cfrq = cfrq+ccr[i]; i=i+1
    i=i+1
    while i < len(ccr) and j < 4 : cyyy = cyyy+ccr[i]; i=i+1; j=j+1
    #
    return [conf, case, cfrq, cyyy]










def mean_3d(XD, LSM, XVOL):
    #
    # XD             : 3D+T array containing data
    # LSM            : 3D land sea mask
    # XVOL           : 3D E1T*E2T*E3T  : 3D mesh sizes
    #
    # RETURN vmean: vector containing mean values for each time

    [ lt, lz, ly, lx ] = nmp.shape(XD)

    if nmp.shape(LSM) != ( lz, ly, lx ):
        print 'ERROR: mean_3d.barakuda_orca.py => XD and LSM do not agree in shape!'
        sys.exit(0)
    if nmp.shape(XVOL) != ( lz, ly, lx ):
        print 'ERROR: mean_3d.barakuda_orca.py => XD and XVOL do not agree in shape!'
        sys.exit(0)

    vmean = nmp.zeros(lt)
    XX = LSM[:,:,:]*XVOL[:,:,:]

    for jt in range(lt):
        vmean[jt] = nmp.sum( XD[jt,:,:,:]*XX ) / nmp.sum( XX )

    return vmean


def mean_2d(XD, LSM, XAREA):
    #
    # XD     : 2D+T array containing data
    # LSM    : 2D land sea mask
    # XAREA  : E1T*E2T, 2D mesh sizes
    #
    # RETURN vmean: the mean value at each record

    [ lt, ly, lx ] = nmp.shape(XD)

    if nmp.shape(LSM) != ( ly, lx ):
        print 'ERROR: mean_3d.barakuda_orca.py => XD and LSM do not agree in shape!'
        sys.exit(0)
    if nmp.shape(XAREA) != ( ly, lx ):
        print 'ERROR: mean_3d.barakuda_orca.py => XD and XAREA do not agree in shape!'
        sys.exit(0)

    vmean = nmp.zeros(lt)
    XX = LSM[:,:]*XAREA[:,:]
    for jt in range(lt):
        vmean[jt] = nmp.sum( XD[jt,:,:]*XX ) / nmp.sum( XX )

    return vmean






def ij_from_xy(xx, yy, xlon, xlat):
    #
    #=============================================================
    # Input:
    #       xx : longitude of searched point (float)
    #       yy : latitude  of searched point (float)
    #       xlon  : 2D array of longitudes of the ORCA grid
    #       xlat  : 2D array of latitudes  of the ORCA grid
    # Output:
    #       [ ji, jj ] i and j indices corresponding
    #=============================================================    
    #
    ji=0 ; jj=0
    #
    if xx < 0.: xx = xx + 360.
    #
    [nj , ni] = xlon.shape
    iwa = nmp.where(xlon < 0.) ; xlon[iwa] = xlon[iwa] + 360. # Want only positive values in longitude:
    #
    # Southernmost latitude of the ORCA domain:
    ji0 = nmp.argmin(xlat[0,:])
    lat_min = xlat[0,ji0]
    ji = ji0
    yy = max( yy, lat_min )
    #print ' lat_min is ', lat_min
    #print ' yy set to ', yy
    # Northernmost latitude of the ORCA domain:    
    ji0 = nmp.argmax(xlat[nj-1,:])
    lat_max = xlat[nj-1,ji0]
    yy = min( yy, lat_max )
    #print " highest latitude:", lat_max
    #
    #im1 = -1000    
    #jm1 = -1000
    #im2 = -1000    
    #jm2 = -1000
    #lfound = False    
    #ji = ji0 ; # Starting with a ji that reaches the northernmost latitude:
    #while not lfound:
    #    #        
    #    # 1/ Finding jj from latitude position and updated ji
    #    jj = nmp.argmin( nmp.abs( xlat[:-2,ji]-yy )) ; # 2 points implied in northfold...
    #    print ' jj =', jj
    #    #
    #    # 2/ Finding ji from longitude position and updated jj:
    #    ji = nmp.argmin( nmp.abs( xlon[jj,:-2]-xx ))
    #    print ' ji =', ji
    #
    #    if (jj-jm1 == 0) and (ji-im1 == 0) : lfound = True
    #    # Prevent infinite loop (oscilating between 2 sets of points):
    #    if (jj-jm2 == 0) and (ji-im2 == 0) : lfound = True
    #    #
    #    im2 = im1
    #    jm2 = jm1
    #    im1 = ji
    #    jm1 = jj
    #
    A = nmp.abs( xlat[:-2,:-2]-yy ) + nmp.abs( xlon[:-2,:-2]-xx )
    (jj,ji) = nmp.unravel_index(A.argmin(), A.shape)
    #
    #print ' Other method: ji, jj =>', ji0, jj0
    #print '  lon, lat =', xlon[jj0,ji0], xlat[jj0,ji0], '\n'
    #
    return ( ji, jj )



def line_vert_hori(x1, x2, y1, y2, xlon, xlat): #, rmin, rmax, dc):
    #
    # LOLO: test and replace coor2ind...
    #
    ji1=0 ; ji2=0 ; jj1=0 ; jj2=0
    lhori = False ; lvert = False
    #
    if y1 == y2: lhori = True
    if x1 == x2: lvert = True
    #
    if not (lhori or lvert) : print 'coor2ind only supports horizontal or vertical sections so far...'; sys.exit(0)
    #
    (ji1,jj1) = ij_from_xy(x1, y1, xlon, xlat)
    (ji2,jj2) = ij_from_xy(x2, y2, xlon, xlat)
    #
    if lhori and (jj1 != jj2): jj2=jj1
    if lvert and (ji1 != ji2): ji2=ji1
    #
    return [ ji1, ji2, jj1, jj2 ]


def coor2ind(x1, x2, y1, y2, xlon, xlat): #, rmin, rmax, dc):
    #
    #=============================================================
    # Input:
    #       x1, x2: longitudes of points 1 and 2 (float)
    #       y1, y2: latitudes  of points 1 and 2 (float)
    #       xlon  : 2D array of longitudes of the ORCA grid
    #       xlat  : 2D array of latitudes  of the ORCA grid
    # Output:
    #       [ ji1, ji2, jj1, jj2 ] i and j indices corresponding
    #=============================================================    
    #
    ji1=0 ; ji2=0 ; jj1=0 ; jj2=0
    lhori = False ; lvert = False
    #
    if x1 < 0.: x1 = x1 + 360.
    if x2 < 0.: x2 = x2 + 360.
    #
    if y1 == y2: lhori = True
    if x1 == x2: lvert = True
    #
    if not (lhori or lvert) : print 'coor2ind only supports horizontal or vertical sections so far...'; sys.exit(0)
    #if lhori: print 'coor2ind horizontal mode not done yet'; sys.exit(0)
    #
    [nj , ni] = xlon.shape
    iwa = nmp.where(xlon < 0.) ; xlon[iwa] = xlon[iwa] + 360. # Want only positive values in longitude:
    #
    jj1 = 0
    # Loop along updated jj1
    for jc in range(5):
        #
        # 1/ Finding ji1 from longitude position
        #    Lookin in regular region of the grid (jj1 pas trop grand, ici = 0):
        dist = 1000.
        for ji in range(ni):
            dd = abs(x1 - xlon[jj1,ji])
            if dd > 360.: dd = dd - 360.
            if dd < dist: dist = dd ; ji1 = ji
        #
        # 2/ Finding jj1 from latitude position and updated ji1
        dist = 1000.
        for jj in range(nj):
            dd = abs(y1 - xlat[jj,ji1])
            if dd < dist: dist = dd ; jj1 = jj
        if y1 == -90.: jj1 = 0
    #
    if lvert:
        ji2 = ji1
        # Neet to find jj2 !
        dist = 1000.
        for jj in range(nj):
            dd = abs(y2 - xlat[jj,ji2])
            if dd < dist: dist = dd ; jj2 = jj
    if y2 ==  90.: jj2 = nj-1
    #
    if lhori:
        jj2 = jj1
        # Neet to find ji2 !
        dist = 1000.
        for ji in range(ni):
            dd = abs(x2 - xlon[jj2,ji])
            if dd > 360.: dd = dd - 360.
            if dd < dist: dist = dd ; ji2 = ji
    #
    return [ ji1, ji2, jj1, jj2 ]





