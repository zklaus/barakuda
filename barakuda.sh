#!/usr/bin/env bash

#==============================================================
#
#                    B A R A K U D A
#
#    An OCEAN MONITORING python environment for NEMO
#
#             L. Brodeau, 2009-2017
#
#===============================================================

export script=barakuda
[ -z ${BARAKUDA_ROOT+x} ] && export BARAKUDA_ROOT=${PWD}

# Display available configs:
list_conf="`\ls ${BARAKUDA_ROOT}/configs/config_*.sh | sed -e "s|${BARAKUDA_ROOT}/configs\/config_||g" -e s/'.sh'/''/g`"
# User configs, potentially in the directory from which barakuda.sh is called:
list_conf+=" `\ls ./config_*.sh 2>/dev/null | sed -e "s|.\/config_||g" -e s/'.sh'/''/g`"


# Important bash functions:
. ${BARAKUDA_ROOT}/src/bash/bash_functions.bash

barakuda_init

while getopts C:R:f:y:c:FeEh option ; do
    case $option in
        C) export CONFIG=${OPTARG} ;;
        R) export EXP=${OPTARG} ;;
        f) export IFREQ_SAV_YEARS=${OPTARG} ;;
        y) export YEAR0=${OPTARG} ; export LFORCE_YINI=true ;;
        c) export EXPREF=${OPTARG} ;;
        F) export LFORCEDIAG=true ;;
        e) export ISTAGE=2 ;;
        E) export ISTAGE=2 ; export l_clim_diag=true ;;
        h)  barakuda_usage; exit ;;
        \?) barakuda_usage; exit ;;
    esac
done

barakuda_check

# Time to source the config file:
if [ -f ./config_${CONFIG}.sh ]; then
    # sourcing local configuration file if present:
    fconfig=./config_${CONFIG}.sh
else
    # sourcing barakuda-distribution configuration file:
    fconfig=${BARAKUDA_ROOT}/configs/config_${CONFIG}.sh
fi
if [ -f ${fconfig} ]; then
    echo "Sourcing configuration file: ${fconfig} !"
    . ${fconfig}
else
    echo "ERROR: cannot find file ${fconfig} !"; exit
fi
echo




# If auto-submit experiment (ece_exp=10) then overides a few functions with:
if [ ${ece_exp} -ge 10 ]; then
    echo "Sourcing ${BARAKUDA_ROOT}/src/bash/bash_functions_autosub.bash"
    . ${BARAKUDA_ROOT}/src/bash/bash_functions_autosub.bash
fi


# If 3D fieds are annual averaged then overides a few functions with:
if [ ! "${ANNUAL_3D}" = "" ]; then
    . ${BARAKUDA_ROOT}/src/bash/bash_functions_1y.bash
fi

# List of CDFTOOLS executables needed for the diagnostics:
export L_EXEC="cdfmaxmoc.x cdfmoc.x cdfvT.x cdftransportiz.x cdficeflux.x cdficediags.x cdfmhst.x cdfsigtrp.x"

barakuda_setup

echo
echo " SETTINGS: "
echo "   *** CONFIG     = ${CONFIG} "
echo "   *** NEMO_OUT_D = ${NEMO_OUT_D} "
echo "   *** CLIM_DIR   = ${CLIM_DIR} "
echo "   *** TMP_DIR    = ${TMP_DIR} "
echo "   *** GRID       = ${ORCA} "
echo "   *** EXP        = ${EXP} "
echo "   *** CPREF      = ${CPREF} "
echo "   *** IFREQ_SAV_YEARS = ${IFREQ_SAV_YEARS} "
echo "   *** NCDF_DIR        = ${NCDF_DIR} "
echo


if [ ${ISTAGE} -eq 1 ]; then
    barakuda_first_last_years ; # look at NEMO files to know what are first and last years available...
    echo ${IFREQ_SAV_YEARS} > ${DIAG_D}/numb_year_per_file.info
    echo ${YEAR_INI}        > ${DIAG_D}/first_year.info
else
    # -> this is stage 2 (plot generation) ISTAGE=2 !
    barakuda_init_plot
fi

cyear_ini=`printf "%04d" ${YEAR_INI}`
cyear_end=`printf "%04d" ${YEAR_END}`

# For proper python executables and scripts to be found:
export PATH=${PYBRKD_EXEC_PATH}:${BARAKUDA_ROOT}/src/bash:${PYTHON_HOME}/bin:${PATH}

#                                   setup over
######################################################################################



jyear=${YEAR_INI}

fcompletion=${DIAG_D}/last_year_done.info
if [ -f ${fcompletion} ]; then jyear=`cat ${fcompletion}`; ((jyear++)); fi

cd ${TMP_DIR}/

barakuda_import_mesh_mask ; # Importing mesh_mask (+basin) files...

if [ ${ISTAGE} -eq 1 ]; then
    # Importing cdftools executables:
    for ex in ${L_EXEC}; do rsync -v ${BARAKUDA_ROOT}/cdftools_light/bin/${ex} . ; done
fi




# L O O P   A L O N G   Y E A R S
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lcontinue=true

if ${LFORCEDIAG}; then lcontinue=false; fi

while ${lcontinue}; do

    export cyear=`printf "%04d" ${jyear}`
    cpf=""
    if [ ${ISTAGE} -eq 1 ] && [ ${ece_exp} -gt 0 ]; then
        iy=$((${jyear}-${YEAR_INI}+1+${YEAR_INI}-${YEAR_INI_F}))
        dir_ece=`printf "%03d" ${iy}`
        echo " *** ${cyear} => dir_ece = ${dir_ece}"
        cpf="${dir_ece}/"
    fi

    i_get_file=0
    if [ $((${jyear}%${IFREQ_SAV_YEARS})) -eq 0 ]; then
        barakuda_check_year_is_complete  ; # lcontinue might be updated to false!
    fi    
    echo; echo
    
    export CRT1M=${CPREF}${cyear}${cmmdd1}_${cyear}${cmmdd2}
    if ${l_y2_j}; then export CRT1M=${CPREF}${cyear}${cmmdd1}_`printf "%04d" $((${jyear}+1))`${cmmdd2} ; fi
    
    if ${lcontinue}; then

        echo "Yeah! Year ${jyear} is saved..."; echo

        
        if [ ${ISTAGE} -eq 2 ]; then
            echo; echo "You cannot create figures and HTML pages yet!"
            echo " => finish treating the results first by launching barakuda.sh without the '-e' switch."
            exit
        fi
        
        echo; echo
        echo "*********************************************************************"
        echo "  Experiment ${EXP}: Will generate diagnostics and data for year ${cyear}..."
        echo "*********************************************************************"
        echo; echo

        barakuda_import_files
        
        # Monthly files to work with for current year:
	# PD: Massive hack: use new CMIP6 files
        ftopa2d1m=${CRT1M}_opa_grid_T_2D.nc
        fuopa2d1m=${CRT1M}_opa_grid_U_2D.nc
        fvopa2d1m=${CRT1M}_opa_grid_V_2D.nc
	ftopa3d1m=${CRT1M}_opa_grid_T_3D.nc
        fuopa3d1m=${CRT1M}_opa_grid_U_3D.nc
        fvopa3d1m=${CRT1M}_opa_grid_V_3D.nc
	ftlim2d1m=${CRT1M}_lim_grid_T_2D.nc
        #fj1m=${CRT1M}_${FILE_ICE_SUFFIX}.nc ; # can be icemod or grid_T ....
        #ff1m=${CRT1M}_${FILE_FLX_SUFFIX}.nc ; # file with surface fluxes
        #
        # Annual files to work with for current year:
        #CRT1Y=`echo ${CRT1M} | sed -e s/"_${TSTAMP}_"/"_${ANNUAL_3D}_"/g`
        #ft1y=${CRT1Y}_grid_T.nc
        #fu1y=${CRT1Y}_grid_U.nc
        #fv1y=${CRT1Y}_grid_V.nc
        #fj1y=${CRT1Y}_${FILE_ICE_SUFFIX}.nc
        #ff1y=${CRT1Y}_${FILE_FLX_SUFFIX}.nc
        CFG3D=${CRT1M}
        #
        # Files that contain the 3D fields (might be monthly or annaual sometimes (when "${ANNUAL_3D}" = "1y")        
        #ft3d=${ft1m}
        #fu3d=${fu1m}
        #fv3d=${fv1m}
        #if [ "${ANNUAL_3D}" = "1y" ]; then
        #    [[ ${NEMO_SAVED_FILES_3D} =~ (^|[[:space:]])"grid_U"($|[[:space:]]) ]] \
        #        && CFG3D=${CRT1Y}; ft3d=${ft1y}; fu3d=${fu1y}; fv3d=${fv1y} \
        #        || echo "...default"
        #    echo ""
        #fi
        fvt=${CFG3D}_VT.nc
        

        # -- time to compute diagnostics --

	# PD: WORKING BUT SLOW AND NOT CLEAR TO WHICH AIM
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # If coupled EC-Earth simu, attempting to compute ocean-averaged fluxes from IFS too (E, P, E-P)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${ece_exp} -eq 2 ] && [ ${NBL} -eq 75 ] && [ ${i_do_ifs_flx} -eq 1 ]; then
            echo; echo; echo "Fluxes of freshwater at the surface from IFS..."
            echo " *** CALLING: extract_ifs_surf_fluxes.sh &"
 	    export PATH=${BARAKUDA_ROOT}/src/bash/misc_ifs:${PATH}
            extract_ifs_surf_fluxes.sh
            pid_flxl=$! ; echo
        fi
	
	# PD: OK!
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Computing time-series of spatially-averaged variables
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_mean} -eq 1 ]; then
            #
            echo; echo "3D-averaging for time series"
            echo " *** CALLING: mean_3d.py ${ftopa3d1m} ${jyear} T &"
            mean_3d.py ${ftopa3d1m} ${jyear} T &
            pid_mn3dt=$! ; echo
            
            echo " *** CALLING: mean_3d.py ${ftopa3d1m} ${jyear} S &"
            mean_3d.py ${ftopa3d1m} ${jyear} S &
            pid_mn3ds=$! ; echo
            
            echo; echo "2D-averaging for time series"
            echo " *** CALLING: mean_2d.py ${ftopa2d1m} ${jyear} &"
            mean_2d.py ${ftopa2d1m} ${jyear} &
            pid_mn2d=$! ; echo

            echo; echo "2D-integration of some surface fluxes over basins"
            echo " *** CALLING: flux_int_basins.py ${ftlim2d1m} ${jyear} &"
            flux_int_basins.py ${ftlim2d1m} ${jyear} &
            pid_ssf=$! ; echo
            #
        fi

	# PD: OK!
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Creating VT file if needed
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_trsp} -gt 0 ] || [ ${i_do_mht} -eq 1 ]; then
            if [ ! -f ${fvt} ]; then
                echo; echo; echo " *** CALLING: ./cdfvT.x ${CFG3D} ${NN_T} ${NN_S} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV} &"
                ./cdfvT.x ${CFG3D} ${NN_T} ${NN_S} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV} &
                pid_vtvt=$! ; echo
            fi
        fi

	# PD: OK!
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2D maps of NEMO - OBS for SST and SSS (for movies)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_movi} -eq 1 ]; then
            echo; echo; echo "2D maps of NEMO - OBS for SST and SSS (for movies)"
            echo " *** CALLING: prepare_movies.py ${ftopa2d1m} ${jyear} sst &"
            prepare_movies.py ${ftopa2d1m} ${jyear} sst &
            pid_movt=$! ; echo
            echo " *** CALLING: prepare_movies.py ${ftopa2d1m} ${jyear} sss &"
            prepare_movies.py ${ftopa2d1m} ${jyear} sss &
            pid_movs=$! ; echo
            echo " *** CALLING: prepare_movies.py ${ftopa2d1m} ${jyear} mld &"
            prepare_movies.py ${ftopa2d1m} ${jyear} mld &
            pid_movm=$! ; echo
            echo " *** CALLING: prepare_movies.py ${ftlim2d1m} ${jyear} ice &"
            prepare_movies.py ${ftlim2d1m} ${jyear} ice &
            pid_movi=$! ; echo
        fi        

	# NOT YET CHECKED
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Computing time-series of spatially-averaged variables
        # on boxes (saving the variable on 2D box too...
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_ssx_box} -eq 1 ]; then
            echo; echo; echo "Box monthly values"
            echo " *** CALLING: ssx_boxes ${ftopa2d1m} ${jyear} ${NN_SST} ${NN_SSS} &"
            ssx_boxes.py ${ftopa2d1m} ${jyear} ${NN_SST} ${NN_SSS} &
            pid_boxa=$! ; echo
        fi

	# PD: OK!
        #~~~~~~~
        #  MOC
        #~~~~~~~
        if [ ${i_do_amoc} -eq 1 ]; then
            echo; echo; echo "MOC"
            rm -f moc.nc *.tmp
            echo " *** CALLING: ./cdfmoc.x ${fvopa3d1m} ${NN_V} ${NN_V_EIV} &"
            ./cdfmoc.x ${fvopa3d1m} ${NN_V} ${NN_V_EIV} &
            pid_amoc=$! ; echo
        fi


	# NOT YET CHECKED
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#  Transport by sigma-class 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_sigt} -eq 1 ]; then
            echo; echo
            if [ ! -f ./dens_section.dat ]; then
                if [ -f ${DENSITY_SECTION_FILE} ]; then
                    echo "Copying ${DENSITY_SECTION_FILE} to here: `pwd` !"; cp ${DENSITY_SECTION_FILE} ./dens_section.dat
                else
                    echo; echo "WARNING: Can't do Transport by sigma-class: ${DENSITY_SECTION_FILE} is missing!!!"
                fi
            fi
            echo " *** CALLING: ./cdfsigtrp.x ${ft3d} ${fu3d} ${fv3d} 24.8 28.6 19 ${jyear} ${DIAG_D}  ${NN_T} ${NN_S} ${NN_U} ${NN_V} &"; echo
            ./cdfsigtrp.x ${ft3d} ${fu3d} ${fv3d} 24.8 28.6 19 ${jyear} ${DIAG_D} ${NN_T} ${NN_S} ${NN_U} ${NN_V}
            pid_sigt=$! ; echo
        fi

        echo
        echo " Gonna wait for level #1 to be done !"
        wait ${pid_vtvt}
        #wait
        echo " .... diag level #1 done...." ; echo        
        echo

        # PD: Extremely slow
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Meridional heat and salt transport
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_mht} -eq 1 ]; then
            echo; echo; echo "Meridional transport of heat and salt"
            fo=${DIAG_D}/merid_transport_T_S_${CONFEXP}.nc
            if [ ! -f ${fvt} ]; then
                echo "WARNING: file ${fvt} is not here, skipping meridional transports section"
            else
                rm -f merid_heat_trp.dat merid_salt_trp.dat
                echo " *** CALLING: ./cdfmhst.x ${fvt} ${fo} ${jyear} &"
                ./cdfmhst.x ${fvt} ${fo} ${jyear} &
                pid_mhst=$! ; echo
            fi
        fi
                
        # PD: OK!
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # VOLUME, HEAT and SALT transports through specified section
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_trsp} -gt 0 ]; then
            echo; echo; echo "Transports of volume, heat and salt through different sections"
            if [ -z ${TRANSPORT_SECTION_FILE} ]; then
                echo "Please specify which TRANSPORT_SECTION_FILE to use into the config file!" ; exit
            fi
            if [ ! -f ./transportiz.dat ]; then
                check_if_file ${TRANSPORT_SECTION_FILE}
                cp ${TRANSPORT_SECTION_FILE} ./transportiz.dat
            fi
            if [ ${i_do_trsp} -eq 1 ]; then z1_trsp="" ; z2_trsp=""; fi
            if [ ! -f ${fvt} ]; then
                echo "WARNING: file ${fvt} is not here, skipping transport section!"
            else
                # Breaking transportiz file in several files so we can go parallel!
                rm -f transportiz_*.dat
                break_downn_section_file.py ./transportiz.dat 2 ; # => Max. 2 sections per file
                list=`\ls transportiz_*.dat`
                for fs in ${list}; do
                    echo " *** CALLING: ./cdftransportiz.x ${fs} ${CFG3D} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV} ${jyear} ${DIAG_D} ${z1_trsp} ${z2_trsp} &"
                    ./cdftransportiz.x ${fs} ${CFG3D} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV} ${jyear} ${DIAG_D} ${z1_trsp} ${z2_trsp} &
                    echo
                done
            fi
        fi

        # NOT WORKING
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Solid freshwater transport through sections due to sea-ice drift
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_ice} -gt 0 ] && [ ${i_do_trsp_ice} -eq 1 ]; then
            echo; echo; echo "Solid freshwater transport through sections due to sea-ice drift"
            if [ -z ${TRANSPORT_SECTION_FILE_ICE} ]; then
                echo "Please specify which TRANSPORT_SECTION_FILE_ICE to use into the config file!" ; exit
            fi
            if [ ! -f ./transport_ice.dat ]; then
                check_if_file ${TRANSPORT_SECTION_FILE_ICE}
                cp ${TRANSPORT_SECTION_FILE_ICE} ./transport_ice.dat
            fi
            #
            # Breaking transportiz file in several files so we can go parallel!
            rm -f transport_ice_*.dat
            break_downn_section_file.py ./transport_ice.dat 2 ; # => Max. 2 sections per file
            list=`\ls transport_ice_*.dat`
            for fs in ${list}; do
                echo " *** CALLING: ./cdficeflux.x ${fs} ${ftlim2d1m} ${NN_ICEF} ${NN_ICEU} ${NN_ICEV} ${NN_ICET} ${jyear} ${DIAG_D} &"
                ./cdficeflux.x ${fs} ${ftlim2d1m} ${NN_ICEF} ${NN_ICEU} ${NN_ICEV} ${NN_ICET} ${jyear} ${DIAG_D} &
                echo
            done
        fi

        # NOT YET CHECKED
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #  Deep Mixed Volume (DMV) on a given box
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ! -z ${i_do_dmv} ] && [ ${i_do_dmv} -gt 0 ]; then

            if [ -z ${FILE_DEF_BOXES} ]; then
                echo "Please specify a FILE_DEF_BOXES to use into the config file!" ; exit
            fi
            echo " *** CALLING: dmv.py ${ft1m} ${cyear} &"
            dmv.py ${ft1m} ${cyear} &
            pid_dmvl=$! ; echo
        fi

	# NOT YET CHECKED
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Budget and other stuffs on a given rectangular box!
        # It provides time-series depending only on time (not depth)
        # budget_rectangle_box.py
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if [ ${i_do_bb} -gt 0 ]; then

            echo; echo; echo "Budget and other stuffs on rectangular boxes!"

            if [ -z ${FILE_DEF_BOXES} ]; then
                echo "Please specify a FILE_DEF_BOXES to use into the config file!" ; exit
            fi
            echo " *** CALLING: budget_rectangle_box.py ${cyear} 100 uv &"
            budget_rectangle_box.py ${cyear} 100 uv &
            pid_bbbb=$! ; echo
        fi

        echo " Gonna wait for level #2 to be done !"
        wait ${pid_amoc} ; # moc needs to be done to call cdfmaxmoc.x ...
        echo
        echo " .... diag level #2 done...." ; echo ; echo


	# PD: OK!
        #~~~~~~~~~~~
        #  Max AMOC
        #~~~~~~~~~~~
        if [ ${i_do_amoc} -eq 1 ]; then
            if [ -z "${LMOCLAT}" ]; then
                echo "AMOC => specify latitude bands with variable LMOCLAT into the config file!!!"; exit
            fi
            for clat in ${LMOCLAT}; do
                cslat=`echo ${clat} | sed -e s/'-'/' '/g`
                echo "  *** ./cdfmaxmoc.x moc.nc atl ${cslat} 500 1500 ${jyear} ${DIAG_D}"
                ./cdfmaxmoc.x moc.nc atl ${cslat} 500 1500 ${jyear} ${DIAG_D} &
            done
        fi

	# PD: OK!
        #~~~~~~~~~
        # SEA-ICE
        #~~~~~~~~~
        if [ ${i_do_ice} -eq 1 ]; then
            echo; echo; echo "Sea-ice extent and volume..." ; rm -f tmp_ice.nc
            echo "ncks  -A -v ${NN_ICEF} ${ftlim2d1m} -o tmp_ice.nc"
            ncks  -A -v ${NN_ICEF} ${ftlim2d1m} -o tmp_ice.nc
            ncrename -v ${NN_ICEF},ice_frac tmp_ice.nc
            coic=""
            if [ -z ${NN_ICET} ]; then
                coic="oic" ; # means only ice concentration available!
            else
                echo "ncks  -A -v ${NN_ICET} ${ftlim2d1m} -o tmp_ice.nc"
                ncks  -A -v ${NN_ICET} ${ftlim2d1m} -o tmp_ice.nc
                ncrename -v ${NN_ICET},ice_thic tmp_ice.nc
            fi
            echo " *** CALLING: ./cdficediags.x tmp_ice.nc ${jyear} ${DIAG_D} ${coic} &"
            ./cdficediags.x tmp_ice.nc ${jyear} ${DIAG_D} ${coic} &

        fi

	# PD: OK!
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Temperature and Salinity on cross meridional/zonal sections
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_sect} -eq 1 ]; then
            check_if_file ${TS_SECTION_FILE}
	    lst_sec=`cat ${TS_SECTION_FILE} | grep -v -E '(^#|EOF)' | awk -F ' ' '{print $1}'`
            echo; echo; echo "Cross-sections on specified transects:"
            echo "${lst_sec}"; echo
            echo " *** CALLING: cross_sections.py ${ftopa3d1m} ${jyear} &"
            cross_sections.py ${ftopa3d1m} ${jyear} &
            echo
        fi

       	# NOT YET CHECKED
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #  Vertical profiles of T,S and density horizontally averaged
        #  on rectangular boxes defined into FILE_DEF_BOXES
        #   => creates time-series function of time and depth
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if [ ${i_do_box_TS_z} -gt 0 ]; then
            if [ -z ${FILE_DEF_BOXES} ]; then
                echo "Please specify a FILE_DEF_BOXES to use into the config file!" ; exit
            fi
            echo " *** CALLING: prof_TS_z_box.py ${cyear} &"
            prof_TS_z_box.py ${cyear} &
            echo;echo
        fi




        # --- end of stuffs that can be launched in bg ---


        echo
        echo " Waiting for backround jobs for current year (${jyear}) !"
        wait  ${pid_mn3dt} ${pid_mn3ds} ${pid_movt} ${pid_movs} ${pid_movm} ${pid_movi} ${pid_mn2d} ${pid_ssf} ${pid_flxl}
        wait
        echo "  Done waiting for year ${cyear} !"
        if [ ${i_do_movi} -eq 1 ]; then rsync -rv movies ${DIAG_D}/ ; fi
        rm -f *.tmp broken_line_* tmp_ice.nc
        rm -f ${CRT1M}_*.nc ${CRT1Y}_*.nc ; #debug
        echo

        echo " ---- DIAGS ARE DONE FOR YEAR ${cyear} ! ---"
        echo "${cyear}" > ${fcompletion}
        echo; echo

    fi ; # if ${lcontinue}; then


    #if ${LFORCE_YEND}; then
    #    if [ ${jyear} -eq ${YEARN} ]; then lcontinue=false; fi
    #fi

    wait

    ((jyear++))


# end loop years...
done ; # while ${lcontinue}; do




# PREPARING HTML PAGE
# ~~~~~~~~~~~~~~~~~~~

l_pclim=false

if [ ${ISTAGE} -eq 2 ]; then

    rm -rf ${DIAG_D}/barakuda/${EXP}
    
    y1=`cat ${DIAG_D}/first_year.info`
    y2=`cat ${DIAG_D}/last_year_done.info`
    nby=$((${y2}-${y1}+1))

    # Agreement between last year from output files and 'fcompletion' file:
    ydum=`cat ${fcompletion}`
    if [ -z ${YEARN} ]; then
        # (if YEARN is not set...)
        if [ ! ${ydum} -eq ${YEAR_END} ]; then
            echo;
            echo "###################################################################"
            echo "ERROR: in ${fcompletion} last_year = ${ydum}"
            echo "         and from stored files files last_year = ${YEAR_END} !"
            echo "###################################################################"
            echo
            exit
        fi
    fi


    

    echo; echo; echo "EXP ${EXP}: creating plots"; echo

    #cd ${BARAKUDA_ROOT}/

    cd ${DIAG_D}/

    echo
    if [ ${i_do_movi} -eq 1 ]; then
        if [ "${iffmpeg_x264}" = "1" ]; then
            # FFMPEG compiled with x264 mp4 support is available on host:
            cd movies/
            for cc in dsst dsss mld icen ices; do
                hh=520 ; # height of image in pixels
                if [ "${cc}" = "icen" ]; then hh=576; fi
                if [ "${cc}" = "ices" ]; then hh=456; fi
                echo " *** CALLING: images2mp4.sh ${cc} ${FIG_FORM} ${hh} 8 18 &"
                images2mp4.sh ${cc} ${FIG_FORM} ${hh} 8 18
                echo
            done
            cd ${DIAG_D}/
        else
            # Faling back on GIF, with 'convert' of imageMagick:
            idelay=$((120-${nby}*8))
            if [ ${idelay} -lt 10 ]; then idelay=10; fi
            rm -f *_${CONFEXP}.gif
            for cc in dsst dsss mld icen ices; do
                echo " *** CALLING: convert -delay ${idelay} -loop 0 movies/${cc}_*.png ${cc}_${CONFEXP}.gif &"
                convert -delay ${idelay} -loop 0 movies/${cc}_*.png ${cc}_${CONFEXP}.gif
                echo
            done
        fi
    fi
    

    # 1D plots to perform
    # ~~~~~~~~~~~~~~~~~~~

    DIAG_1D_LIST=""

    if [ ${i_do_mean} -eq 1 ]; then
        DIAG_1D_LIST="${DIAG_1D_LIST} 3d_so mean_sos 3d_thetao mean_tos mean_zos"
        if [ ! "${NN_MLD}"  = "X" ]; then DIAG_1D_LIST="${DIAG_1D_LIST} mean_mld"; fi
        if [ ! "${NN_QNET}" = "X" ]; then DIAG_1D_LIST="${DIAG_1D_LIST} mean_htf"; fi
        if [ ! "${NN_FWF}"  = "X" ]; then DIAG_1D_LIST="${DIAG_1D_LIST} mean_fwf"; fi
    fi
    if [ ${i_do_amoc} -eq 1 ]; then DIAG_1D_LIST="${DIAG_1D_LIST} amoc";        fi
    if [ ${i_do_ice}  -eq 1 ]; then DIAG_1D_LIST="${DIAG_1D_LIST} seaice";      fi
    if [ ${i_do_trsp} -gt 0 ] || [ ${i_do_trsp_ice} -eq 1 ]; then DIAG_1D_LIST="${DIAG_1D_LIST} transport_sections" ; fi

    dy=$((${YEAR_END}-${YEAR_INI}+1)) ; export YF2=$((${YEAR_END}+1))


    # Doing 1D plots
    # ~~~~~~~~~~~~~~    

    echo ; echo; echo "Going to perform the following 1D plots:"
    echo "    => ${DIAG_1D_LIST}"; echo

    for fd in ${DIAG_1D_LIST}; do
        echo " *** CALLING: plot_time_series.py ${fd}"
        plot_time_series.py ${fd} ; echo
        echo
    done
    echo ; echo ; echo

    if [ ${i_do_mean} -eq 1 ]; then
        
        # 5-month-running mean SST anomaly over Nino region 3.4 graph:
        echo " *** CALLING: plot_enso.py Nino34_${CONFEXP}.nc ${NN_SST}"
        plot_enso.py Nino34_${CONFEXP}.nc ${NN_SST}
        echo; echo

        # Hovmuller of temperature and salinity
        echo " *** CALLING: plot_hovm_tz.py"
        plot_hovm_tz.py
        echo; echo
        
        if [ ${nby} -ge 70 ]; then
            # AMO aka 11-year-running mean SST anomaly over North Atlantic (0-70N)
            echo " *** CALLING: plot_amo.py mean_SST_NAtl_${CONFEXP}.nc ${NN_SST}"
            plot_amo.py mean_SST_NAtl_${CONFEXP}.nc ${NN_SST}
            echo; echo
        fi        
    fi


    if [ ${i_do_sigt} -eq 1 ]; then
        # Transport by sigma-class
        echo " *** CALLING: plot_trsp_sigma.py"
        plot_trsp_sigma.py
        echo; echo; echo
    fi


    if [ ${i_do_mht} -eq 1 ]; then
        #
        # Hovmullers of advective meridional heat/salt transport
        echo; echo
        echo " *** CALLING: plot_hovm_merid_trsp.py"
        plot_hovm_merid_trsp.py
        echo; echo; echo
        #
    fi



    echo; echo; echo




    if ${l_clim_diag} ; then

        ###########################################################################
        # Climatology over X years (12-month file average of X consecutive years)
        #   => has to be built with the 'build_clim.sh' script
        ###########################################################################

        echo; echo; echo "Checking for presence of ${DIAG_D}/clim/last_clim..."
        if [ -f ${DIAG_D}/clim/last_clim ]; then
            cat ${DIAG_D}/clim/last_clim
            export CLIM_PER=`cat ${DIAG_D}/clim/last_clim`
            ftcli=${DIAG_D}/clim/mclim_${CONFEXP}_${CLIM_PER}_grid_T.nc4
            ficli=${DIAG_D}/clim/mclim_${CONFEXP}_${CLIM_PER}_${FILE_ICE_SUFFIX}.nc4
            fcsbc=${DIAG_D}/clim/mclim_${CONFEXP}_${CLIM_PER}_${FILE_FLX_SUFFIX}.nc4
            fclvt=${DIAG_D}/clim/aclim_${CONFEXP}_${CLIM_PER}_VT.nc4
            fcmoc=${DIAG_D}/clim/aclim_${CONFEXP}_${CLIM_PER}_MOC.nc4
            fcpsi=${DIAG_D}/clim/aclim_${CONFEXP}_${CLIM_PER}_PSI.nc4
            fccrl=${DIAG_D}/clim/aclim_${CONFEXP}_${CLIM_PER}_TCURL.nc4
            iclyear=`echo ${CLIM_PER} | sed -e s/'-'/' '/g`
        else
            echo; echo "ERROR => you set l_clim_diag to true (probably by using '-E') but no file 'last_clim' was found in:"
            echo "               ${DIAG_D}/clim/"
            echo "               => then use '-e' instead!"; echo
            exit
        fi

        echo; echo; echo "Checking for presence of ${ftcli}..."
        if [ -f ${ftcli} ]; then
            echo; echo;
            echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo "*   Climatologies found !!!"
            echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            echo
            echo "  => for years ${CLIM_PER}" ; echo "  => using ${ftcli}"

            list_comp_2d="OBS"
            l_pclim=true
            lcomp_to_exp=false

            if [ ! -z ${EXPREF} ]; then
                lcomp_to_exp=true
                list_comp_2d="OBS ${EXPREF}"
                # Must check if climatology for exp ${EXPREF} is there:
                fclim_ref=`echo "${ftcli}" | sed -e "s|${EXP}|${EXPREF}|g"`
                check_if_file ${fclim_ref}
                echo "Going to compare also against exp ${fclim_ref}!"
                echo
            fi

            echo; echo
            for ff in ${F_T_OBS_3D_12} ${F_S_OBS_3D_12} ${F_SST_OBS_12}; do check_if_file ${ff} "name:${ff}"; done
            if [ ${i_do_ice} -gt 0 ]; then check_if_file ${F_ICE_OBS_12}    "name:${F_ICE_OBS_12}" ; fi
            echo; echo


            #######################################
            # Diags that don't imply a comparison #
            #######################################

            export COMP2D="OBS"
            
            # Lat-Depth AMOC
            # ~~~~~~~~~~~~~~
            if [ -f ${fcmoc} ]; then
                echo; echo
                echo " Ploting lat-depth MOC !"
                cd ${DIAG_D}/
                export DIRS_2_EXP="${DIRS_2_EXP} moc"
                rm -rf moc; mkdir moc; cd moc/
                echo; echo; echo " *** CALLING: moc.py ${iclyear}"
                moc.py ${iclyear}
                cd ../
                echo
            fi

            # March Mixed layer depth in Nordic Seas
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if [ `ipresent_var_in_ncf ${ftcli} ${NN_MLD}` -eq 1 ]; then
                echo; echo
                echo " Performing 2D mapping of March Mixed layer depth in Nordic Seas"
                cd ${DIAG_D}/
                export DIRS_2_EXP="${DIRS_2_EXP} mld"
                rm -rf mld; mkdir mld; cd mld/
                echo; echo; echo " *** CALLING: mld.py ${iclyear}"; echo
                mld.py ${iclyear}
                cd ../
                echo
            fi

            # Sea-ice extent stereographic polar projection South and North
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if [ ${i_do_ice}  -gt 0 ] && [ `ipresent_var_in_ncf ${ficli} ${NN_ICEF}` -eq 1 ]; then
                echo; echo
                echo " Performing 2D Sea-ice extent stereographic polar projection South and North"
                cd ${DIAG_D}/
                export DIRS_2_EXP="${DIRS_2_EXP} sea_ice"
                rm -rf sea_ice; mkdir sea_ice; cd sea_ice/
                echo; echo; echo " *** CALLING: ice.py ${iclyear}"; echo
                ice.py ${iclyear}
                cd ../
                echo
            fi

            # Sea-surface height
            # ~~~~~~~~~~~~~~~~~~
            if [ `ipresent_var_in_ncf ${ftcli} ${NN_SSH}` -eq 1 ]; then
                echo; echo; echo " SSH map"
                export DIRS_2_EXP="${DIRS_2_EXP} ssh"
                cd ${DIAG_D}/
                rm -rf ssh; mkdir ssh; cd ssh/
                echo " *** CALLING: ssh.py ${iclyear}"; echo
                ssh.py ${iclyear} 
                cd ../ ;  echo
            else
                echo; echo "WARNING: did not find ${NN_SSH} into ${ftcli} !!!!"; echo
            fi
            
            # Wind-stress module and curl
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if [ -f ${fccrl} ]; then
                echo; echo; echo " Wind-stress module and curl maps"
                if [ -x ${NN_TAUM} ]; then echo "ERROR: define variable NN_TAUM in ${fconfig}! ('X' if not present in NEMO output)"; exit; fi
                export DIRS_2_EXP="${DIRS_2_EXP} wind"
                cd ${DIAG_D}/
                rm -rf wind; mkdir wind; cd wind/
                echo " *** CALLING: wind.py ${iclyear}"; echo
                wind.py ${iclyear}
                cd ../ ;  echo
            else
                echo
                echo "WARNING: did not find file ${fccrl} !!!"
                echo "         or did not find ${NN_TAUM} into ${fcsbc} !!!!"
                echo
            fi
            

            # Zonally-averaged surface heat fluxes
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ lolo
            if [ -f ${fcsbc} ]; then
                echo; echo; echo " Zonally-averaged surface heat fluxes"
                if [ -x ${NN_QSOL} ]; then echo "ERROR: define variable NN_QSOL in ${fconfig}! ('X' if not present in NEMO output)"; exit; fi
                export DIRS_2_EXP="${DIRS_2_EXP} sfluxes"
                cd ${DIAG_D}/
                rm -rf sfluxes; mkdir sfluxes; cd sfluxes/
                echo " *** CALLING: sfluxes.py ${iclyear}"; echo
                sfluxes.py ${iclyear}
                cd ../ ;  echo
            else
                echo
                echo "WARNING: did not find file ${fcsbc} !!!"
                echo "         or did not find ${NN_QSOL} into ${fcsbc} !!!!"
                echo
            fi
            

            ##################################################
            # Diags that imply a comparison against "COMP2D" #
            ##################################################

            cd ${DIAG_D}/
            rm -rf temp_sal

            for COMP2D in ${list_comp_2d}; do

                export COMP2D=${COMP2D}
                echo; echo; echo "Clim. comparisons against ${COMP2D}"
               
                if [ "${COMP2D}" = "${EXPREF}" ]; then
                    export F_T_OBS_3D_12=${fclim_ref}; check_if_file ${F_T_OBS_3D_12} "name:F_T_OBS_3D_12"
                    export F_S_OBS_3D_12=${fclim_ref}; check_if_file ${F_S_OBS_3D_12} "name:F_S_OBS_3D_12"
                    export F_SST_OBS_12=${fclim_ref} ; check_if_file ${F_SST_OBS_12}  "name:F_SST_OBS_12"
                    if [ ${i_do_ice}  -gt 0 ]; then export F_ICE_OBS_12=${fclim_ref}   ; check_if_file ${F_ICE_OBS_12}    "name:F_ICE_OBS_12"; fi
                fi
                
                # Temperature and Salinity
                # ~~~~~~~~~~~~~~~~~~~~~~~~
                echo; echo
                echo " Creating maps and cross-sections (i_do_sect) of Temperature and Salinity"
                cd ${DIAG_D}/
                export DIRS_2_EXP="${DIRS_2_EXP} temp_sal"
                DIRS_2_EXP_RREF="${DIRS_2_EXP_RREF} temp_sal"
                mkdir -p temp_sal; cd temp_sal/
                echo; echo; echo " *** CALLING: temp_sal.py ${iclyear} &"; echo
                temp_sal.py ${iclyear}
                cd ../
                echo
                
            done ; # for COMP2D in ${list_comp_2d}; do

            wait
            
        else
            echo; echo
            echo " No Climatologies found ...";
            echo "   => you can use 'build_clim.sh' to build a climato of your experiment"
            echo; echo
        fi

    fi ; # if ${l_clim_diag}


    
    wait


    # Time for HTML stuff!

    export HTML_DIR=${DIAG_D}/${EXP}/barakuda
    mkdir -p ${HTML_DIR}
    
    cd ${DIAG_D}/

    # Moving all figures to HTML_DIR:
    for fp in ${FIG_FORM} svg mp4 gif; do mv -f *.${fp} ${HTML_DIR}/ >/dev/null 2>/dev/null ; done
    mv -f ./merid_transport/*.${FIG_FORM} ${HTML_DIR}/ >/dev/null 2>/dev/null
    mv -f ${DIAG_D}/movies/movie_*.mp4    ${HTML_DIR}/ >/dev/null 2>/dev/null
    
    # Sign with logo?
    #cd ${HTML_DIR}/
    #list=`\ls *.png`
    #for ff in ${list}; do
    #    sign_image ${ff} n_${ff}
    #    mv -f n_${ff} ${ff}
    #done
    
    cd ${DIAG_D}/

    . ${BARAKUDA_ROOT}/src/bash/build_html.bash

    wait

    # Building main index.html 
    build_index_html
    
    # If climatology built, sub 2D html pages
    if ${l_pclim}; then
        build_sub_html
    fi
    
    #==================================================================
    

    wait ; # likely waiting for the creation of the GIFs....

    echo; echo

    cp ${BARAKUDA_ROOT}/src/html/conf_*.html ${HTML_DIR}/
    
    if [ ${ece_exp} -eq 0 ]; then
        cp -L ${BARAKUDA_ROOT}/src/html/logo.*g      ${HTML_DIR}/
    else
        cp -L ${BARAKUDA_ROOT}/src/html/logo_ece.svg ${HTML_DIR}/logo.svg
        cp -L ${BARAKUDA_ROOT}/src/html/logo_ece.png ${HTML_DIR}/logo.png
    fi
    
    mv -f index.html namelist*.html ${HTML_DIR}/

    cp -r ${DIRS_2_EXP} ${HTML_DIR}/ >/dev/null 2>/dev/null

    echo; echo; echo

    if [ ${ihttp} -eq 1 ]; then
        echo "Preparing to export to remote host!"; echo
        #send_dir=`basename ${HTML_DIR}`
        #tar cvf ${send_dir}.tar ${send_dir}
        #ssh ${RUSER}@${RHOST} "mkdir -p ${RWWWD}"
        #scp ${send_dir}.tar ${RUSER}@${RHOST}:${RWWWD}/
        #ssh ${RUSER}@${RHOST} "cd ${RWWWD}/; rm -rf ${send_dir}; tar xf ${send_dir}.tar 2>/dev/null; rm -f ${send_dir}.tar; \
        #                       chmod -R a+r ${send_dir}; cd ${send_dir}/"
        #echo; echo
        #echo "Diagnostic page installed on  http://${RHOST}${RWWWD}/${send_dir}/ !"
        #echo "( Also browsable on local host in ${DIAG_D}/${send_dir} )"
        #rm -f ${send_dir}.tar

	#PD modify to use ectrans on ECMWF
	tar cvfz barakuda-${EXP}.tar.gz ${EXP}/
	ectrans -remote $RHOST -source barakuda-${EXP}.tar.gz -verbose

    else
        if [ ${ihttp} -eq 0 ]; then
            echo "Diagnostic page installed in ${HTML_DIR}/"
            echo " => you can view it with a web browser..."
        else
            echo "Error: \"ihttp\" must be either 0 or 1 !"
        fi
    fi

    echo; echo

    rm -rf *.eps


else
    echo
    echo "Diagnostics built and saved! (${DIAG_D})"
    echo "Re-Run \"${script}.sh\" adding the  \"-e\" or \"-E\" switch to create figure and HTML page..."
    echo
fi

rm -rf ${TMP_DIR} 2>/dev/null ; #debug

echo
