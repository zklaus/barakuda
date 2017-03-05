#!/usr/bin/env python

#       B a r a K u d a
#
#     Compare time-series between some experiments!
#
#       L. Brodeau, 2013

import sys
import os
import numpy as nmp

from netCDF4 import Dataset

import barakuda_ncio as bn
import barakuda_orca as bo
import barakuda_plot as bp
import barakuda_tool as bt

DEFAULT_LEGEND_LOC = 'out'

iamoc  = 1
i2dfl  = 1
i3dfl  = 1
imld   = 1
iice   = 1
itrsp  = 1
ifwf   = 1  ; # freshwater fluxes at the surface


venv_needed = {'LIST_EXPS','DIAG_DIR','CONF','FIG_FORMAT', 'BM_FILE', \
               'NN_SST','NN_SST','NN_SSS','NN_SSH','NN_T','NN_S','NN_MLD', \
               'TRANSPORT_SECTION_FILE','LMOCLAT','i_do_ifs_flx'}

vdic = bt.check_env_var(sys.argv[0], venv_needed)

cd_diag = vdic['DIAG_DIR']
cffig   = vdic['FIG_FORMAT']


i_do_ifs_flx = int(vdic['i_do_ifs_flx'])

narg = len(sys.argv)
if narg != 3: print 'Usage: '+sys.argv[0]+' <first_year> <last_year>'; sys.exit(0)
cy1 = sys.argv[1] ; y1 = int(cy1)
cy2 = sys.argv[2] ; y2 = int(cy2)

nb_years = y2 - y1 + 1


clist_exps = vdic['LIST_EXPS'].split()
clist_confexps = []

for cexp in clist_exps:
    clist_confexps.append(vdic['CONF']+'-'+cexp)

print sys.argv[0]+': will compare following experiments: '; print clist_confexps
print ' ... saved into '+cd_diag+'\n'

nbexp = len(clist_confexps)



ittic = bt.iaxe_tick(nb_years)

Vt = nmp.zeros(        nb_years)
Xf = nmp.zeros((nbexp, nb_years))


def test_nb_mnth_rec(nbmn, nbyr, cnd):
    print ' *** nb. mnth. records =', nbmn
    if nbmn%12 != 0:
        print 'ERROR: compare_time_series.py => number of monthly records is not a multile of 12 in the netcdf file! diag = '+cnd
        sys.exit(0)
    if nbmn/12 > nbyr:
        print 'ERROR: compare_time_series.py => too many monthly records in netcdf file! diag = '+cnd
        print '     number of expected monthy records =', nbyr*12
        print '                          number found =', nbmn
        sys.exit(0)
    return






# Only one column to read:
##########################

if i2dfl == 1 or i3dfl == 1 or imld == 1 :
    list_basin_names, list_basin_lgnms = bo.get_basin_info(vdic['BM_FILE'])

if i2dfl == 1:

    vvar  = [ vdic['NN_SSH'], vdic['NN_SSS'], vdic['NN_SST'] ]
    vname = [ 'SSH'     ,  'SSS'      , 'SST'     ]
    vunit = [ r'm'    ,  r'PSU'   , r'$^{\circ}$C']


    jvar=0
    for cvar in vvar:
        cdiag = 'mean_'+cvar
        print '\n Treating '+cdiag
        
        joce = 0
        for cocean in list_basin_names :
            Xf[:,:] = 0. ; jexp = 0
            for confexp in clist_confexps:

                cf_in = cd_diag+'/'+confexp+'/'+cdiag+'_'+confexp+'_'+cocean+'.nc'
                vt0, vd0 = bn.read_1d_series(cf_in, cvar, cv_t='time', l_return_time=True)
                nbm = len(vt0)
                test_nb_mnth_rec(nbm, nb_years, cdiag)
                VY, FY = bt.monthly_2_annual(vt0, vd0)
                Vt[:nbm/12]   = VY[:]
                Xf[jexp,:nbm/12] = FY[:] ; Xf[jexp,nbm/12:] = -999.
                jexp = jexp + 1

            bp.plot("1d_multi")(Vt[:], Xf[:,:], clist_exps, cfig_type=cffig,
                                cfignm=cdiag+'_comparison_'+cocean, dt=ittic,
                                loc_legend=DEFAULT_LEGEND_LOC, cyunit=vunit[jvar],
                                ctitle = vname[jvar]+', '+list_basin_lgnms[joce], ymin=0, ymax=0)
            joce = joce + 1
        jvar = jvar+1











if imld == 1:

    cvar  = vdic['NN_MLD']
    cdiag = 'mean_mld'
    vname = [ r'Mixed layer depth' ]
    lplot = True

    print '\n Treating '+cdiag
    Xf[:,:] = 0
    
    joce = 0
    for coce in list_basin_names:
        jexp = 0
        for confexp in clist_confexps:
            cf_in = cd_diag+'/'+confexp+'/mean_'+cvar+'_'+confexp+'_'+coce+'.nc'
            bt.chck4f(cf_in, script_name='compare_time_series.py')
            vt0, vd0 = bn.read_1d_series(cf_in, cvar, cv_t='time', l_return_time=True)
            nbm = len(vt0)
            test_nb_mnth_rec(nbm, nb_years, cdiag)
            VY, FY = bt.monthly_2_annual(vt0, vd0)
            Vt[:nbm/12]   = VY[:]
            Xf[jexp,:nbm/12] = FY[:] ; Xf[jexp,nbm/12:] = -999.
            lplot = lplot and lplot
            jexp = jexp + 1

        if lplot:
            bp.plot("1d_multi")(Vt[:], Xf[:,:], clist_exps, cfig_type=cffig,
                                cfignm=cdiag+'_'+coce+'_comparison', dt=ittic,
                                loc_legend=DEFAULT_LEGEND_LOC, cyunit='m',
                                ctitle = 'Mixed layer depth, '+list_basin_lgnms[joce],
                                ymin=0, ymax=0)
        joce = joce+1











# Several columns to read
#=========================

if i3dfl == 1:

    vvar  = [ vdic['NN_S'],     vdic['NN_T'] ]
    vname = [ 'Salinity' , 'Potential Temperature' ]
    vunit = [ r'PSU'     ,  r'$^{\circ}$C' ]

    vdepth_infil = [ '0-bottom', '0-100', '100-1000', '1000-bottom' ] ; # for 3d_thetao and 3d_so
    vdepth_range = [ 'All', '0m-100m', '100-1000m', '1000m-bottom' ] ; # for 3d_thetao and 3d_so

    jvar=0
    for cvar in vvar:
        cdiag = '3d_'+cvar
        print '\n Treating '+cdiag

        joce = 0
        for cocean in list_basin_names :

            # along the 4 columns of temperature
            idepth=0
            for cdepth in vdepth_range:
                cdif = vdepth_infil[idepth]
                Xf[:,:] = 0.
                jexp = 0
                for confexp in clist_confexps:

                    cf_in = cd_diag+'/'+confexp+'/'+cdiag+'_'+confexp+'_'+cocean+'.nc'
                    vt0, vd0 = bn.read_1d_series(cf_in, cvar+'_'+cdif, cv_t='time', l_return_time=True)
                    nbm = len(vt0)
                    test_nb_mnth_rec(nbm, nb_years, cdiag)
                    VY, FY = bt.monthly_2_annual(vt0, vd0)
                    Vt[:nbm/12]   = VY[:]
                    Xf[jexp,:nbm/12] = FY[:]  ; Xf[jexp,nbm/12:] = -999.
                    jexp = jexp + 1

                bp.plot("1d_multi")(Vt[:], Xf[:,:], clist_exps, cfig_type=cffig,
                                    cfignm=cdiag+'_comparison_'+cocean+'_'+cdepth, dt=ittic, loc_legend=DEFAULT_LEGEND_LOC,
                                    cyunit=vunit[jvar], ctitle = vname[jvar]+', '+list_basin_lgnms[joce]+', depth range = '+cdepth, ymin=0, ymax=0)

                idepth = idepth + 1
                
            joce = joce + 1

        jvar = jvar+1






# Sea-ice
#########

if iice == 1:

    vvar  = [ 'volu'           ,     'area'       ]
    vvnm  = [ 'Sea-ice Volume' , 'Sea-ice Extent' ]
    vunit = [ r'$10^3km^3$'    , r'$10^6$km$^2$'  ]

    vpole = [ 'Arctic', 'Antarctic' ]
    vlab  = [ 'ne'    , 'se'        ]

    jvar = -1
    for cvar in vvar:
        jvar = jvar + 1

        vmnth = [            2          ,              8             ]
        vname = [ vvnm[jvar]+' in March', vvnm[jvar]+' in September' ]

        for jdiag in range(len(vmnth)):
            print '\n Treating '+vname[jdiag]

            # along the 2 columns of Arcit/Antarctic
            ipole = 0
            for cpole in vpole:
                Xf[:,:] = 0.
                jexp = 0
                for confexp in clist_confexps:
                    cf_in = cd_diag+'/'+confexp+'/seaice_diags.nc'
                    vt0, vd0 = bn.read_1d_series(cf_in, cvar+'_'+vlab[ipole], cv_t='time', l_return_time=True)
                    nbm = len(vt0)
                    test_nb_mnth_rec(nbm, nb_years, cdiag)
                    Vt[:nbm/12]   = vt0[vmnth[jdiag]::12]
                    Xf[jexp,:nbm/12] = vd0[vmnth[jdiag]::12] ; Xf[jexp,nbm/12:] = -999.
                    jexp = jexp + 1

                cdiag = 'seaice_'+cvar
                cmnth = '%2.2i'%(vmnth[jdiag]+1)
                bp.plot("1d_multi")(Vt, Xf, clist_exps, cfig_type=cffig,
                                    cfignm=cdiag+'_m'+str(cmnth)+'_comparison_'+cpole, dt=ittic, loc_legend=DEFAULT_LEGEND_LOC,
                                    cyunit=vunit[jdiag], ctitle = vname[jdiag]+', '+cpole, ymin=0, ymax=0)

                ipole = ipole + 1










# Transport through sections
############################

if itrsp == 1:

    print '\nUsing TRANSPORT_SECTION_FILE = '+vdic['TRANSPORT_SECTION_FILE']
    list_sections = bt.get_sections_from_file(vdic['TRANSPORT_SECTION_FILE'])
    print 'List of sections to treat: ', list_sections
    nbsect = len(list_sections)

    vstuff = [ 'volume', 'heat' , 'salt' ]
    vunit  = [ 'Sv'    , 'PW'   , 'kt/s' ]


    jexp = 0
    for confexp in clist_confexps:

        jsect=0
        for csect in list_sections:
            print '\n Treating transports through '+csect

            cf_in = cd_diag+'/'+confexp+'/transport_sect_'+csect+'.nc' ; bt.chck4f(cf_in, script_name='compare_time_series.py')
            id_in = Dataset(cf_in)
            if jsect == 0:
                if jexp == 0:
                    vyear = nmp.zeros(nb_years)
                    Xtrsp = nmp.zeros((nbexp,nbsect,3,nb_years))

                Vt_t = id_in.variables['time'][:] ; nbm = len(Vt_t) ; nby = nbm/12
                test_nb_mnth_rec(nbm, nb_years, cdiag)

            Xtrsp[jexp,jsect,:,:] = -999.
            vyear[:nby], Xtrsp[jexp,jsect,0,:nby] = bt.monthly_2_annual(Vt_t, id_in.variables['trsp_volu'][:nbm])
            vyear[:nby], Xtrsp[jexp,jsect,1,:nby] = bt.monthly_2_annual(Vt_t, id_in.variables['trsp_heat'][:nbm])
            vyear[:nby], Xtrsp[jexp,jsect,2,:nby] = bt.monthly_2_annual(Vt_t, id_in.variables['trsp_salt'][:nbm])

            id_in.close()

            jsect = jsect + 1
        jexp = jexp + 1

    # All data read!

    jsect=0
    for csect in list_sections:
        jstuff = 0
        for cstuff in vstuff:

            bp.plot("1d_multi")(vyear[:], Xtrsp[:,jsect,jstuff,:], clist_exps, cfig_type=cffig,
                                cfignm='transport_'+cstuff+'_'+csect+'_comparison', dt=ittic, loc_legend=DEFAULT_LEGEND_LOC,
                                cyunit=vunit[jstuff], ctitle = 'Transport of '+cstuff+' through section '+csect,
                                ymin=0, ymax=0)

            jstuff = jstuff + 1
        jsect = jsect+1







# AMOC
if iamoc == 1:

    list_lat = vdic['LMOCLAT'].split() ; nblat = len(list_lat)
    print '\n AMOC: '+str(nblat)+' latitude bands!'

    nbm_prev = 0
    jexp = 0
    for confexp in clist_confexps:

        jl = 0
        for clr in list_lat:
            [ c1, c2 ] = clr.split('-') ; clat_info = '+'+c1+'N+'+c2+'N'
            cf_in = cd_diag+'/'+confexp+'/max_moc_atl_'+clat_info+'.nc' ; bt.chck4f(cf_in, script_name='compare_time_series.py')

            id_in = Dataset(cf_in)
            if jl == 0:
                if jexp == 0:
                    vyear = nmp.zeros(nb_years)
                    Xamoc = nmp.zeros((nbexp, nblat, nb_years))
                Vt_t = id_in.variables['time'][:] ; nbm = len(Vt_t) ; nby = nbm/12
                test_nb_mnth_rec(nbm, nb_years, cdiag)

            Xamoc[jexp,jl,:] = -999.
            vyear[:nby], Xamoc[jexp,jl,:nby] = bt.monthly_2_annual(Vt_t[:nbm], id_in.variables['moc_atl'][:nbm])
            id_in.close()

            jl = jl + 1
        jexp = jexp + 1


    jl = 0
    for clr in list_lat:

        bp.plot("1d_multi")(vyear[:], Xamoc[:,jl,:], clist_exps, cfig_type=cffig,
                            cfignm='AMOC_'+clr+'_comparison', loc_legend=DEFAULT_LEGEND_LOC,
                            dt=ittic, cyunit='Sv', ctitle = 'AMOC ('+clr+')', ymin=0, ymax=0)

        jl = jl + 1






# Freshwater Fluxes
if ifwf == 1:

    vvar  = [ 'EmPmR', 'E',           'R'      ,   'P'   , 'ICalv'  ]
    vname = [ 'E-P-R', 'Evaporation', 'Runoffs', 'Precip', 'Calving']
    vunit = [ r'Sv'  , r'Sv',        r'Sv'     ,  r'Sv'  ,  r'Sv'   ]

    jvar=0
    for cvar in vvar:
        cdiag = cvar
        print '\n Treating FWF : '+cdiag

        # Testing only on the first experiment if cvar has been written:
        id_test = Dataset(cd_diag+'/'+clist_confexps[0]+'/mean_fwf_'+clist_confexps[0]+'_GLO.nc')
        list_variables = id_test.variables.keys()
        
        if cvar in list_variables:
            jexp=0
            for confexp in clist_confexps:
    
                cf_in = cd_diag+'/'+confexp+'/mean_fwf_'+confexp+'_GLO.nc'
    
                vt0, vd0 = bn.read_1d_series(cf_in, cvar, cv_t='time', l_return_time=True)
                nbm = len(vt0)
                test_nb_mnth_rec(nbm, nb_years, cdiag)
    
                VY, FY = bt.monthly_2_annual(vt0, vd0)
                Vt[:nbm/12]      = VY[:]
                Xf[jexp,:nbm/12] = FY[:] ; Xf[jexp,nbm/12:] = -999.
                jexp = jexp + 1
    
            bp.plot("1d_multi")(Vt[:], Xf[:,:], clist_exps, cfig_type=cffig,
                                cfignm='FWF_'+cdiag+'_comparison', dt=ittic, loc_legend=DEFAULT_LEGEND_LOC,
                                cyunit=vunit[jvar], ctitle = vname[jvar]+' flux integrated over oceans (NEMO)',
                                ymin=0, ymax=0)
        jvar = jvar+1



    if i_do_ifs_flx == 1:

        # Checking if there are filef for IFS:
        l_fwf_ifs = True  ;  jexp=0
        for confexp in clist_confexps:
            cf_in = cd_diag+'/'+confexp+'/mean_fwf_IFS_'+clist_exps[jexp]+'_GLO.nc'
            print '  *** Checking for the existence of '+cf_in
            if os.path.exists(cf_in):
                print "  *** IFS FWF files found!"
                l_fwf_ifs = True and l_fwf_ifs
            jexp=jexp+1
    
    
        if l_fwf_ifs:
    
            co = ' oceans (IFS)'
            cl = ' land (IFS)'
            vvar  = [ 'flx_e_sv', 'flx_p_sv'  , 'flx_emp_sv', 'flx_e_land_sv', 'flx_p_land_sv'  , 'flx_emp_land_sv' ]
            vname = [ 'E'+co   , 'Precip'+co, 'E-P'+co    , 'E'+cl        , 'Precip'+cl     , 'E-P'+cl          ]
            vunit = [ r'Sv'     ,  r'Sv'      ,  r'Sv'      , r'Sv'          ,  r'Sv'           ,  r'Sv'            ]
            vstit = [ co, co, co, cl, cl, cl ]
    
            jvar=0
            for cvar in vvar:
                cdiag = cvar
                print '\n Treating IFS FWF : '+cdiag
    
                jexp=0
                for confexp in clist_confexps:
                    cf_in = cd_diag+'/'+confexp+'/mean_fwf_IFS_'+clist_exps[jexp]+'_GLO.nc'
    
                    vt0, vd0 = bn.read_1d_series(cf_in, cvar, cv_t='time', l_return_time=True)
                    nbm = len(vt0)
                    test_nb_mnth_rec(nbm, nb_years, cdiag)
    
                    VY, FY = bt.monthly_2_annual(vt0, vd0)
                    Vt[:nbm/12]      = VY[:]
                    Xf[jexp,:nbm/12] = FY[:]  ; Xf[jexp,nbm/12:] = -999.
                    jexp = jexp + 1
    
                bp.plot("1d_multi")(Vt[:], Xf[:,:], clist_exps, cfig_type=cffig,
                                    cfignm='FWF_'+cdiag+'_IFS_comparison', dt=ittic, loc_legend=DEFAULT_LEGEND_LOC,
                                    cyunit=vunit[jvar], ctitle = vname[jvar]+' flux integrated over'+vstit[jvar], ymin=0, ymax=0)
    
                jvar = jvar+1




print  '\n\n'+sys.argv[0]+' done...\n'
