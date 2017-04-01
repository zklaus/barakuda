#!/bin/bash

#lolo_tmp
#export i_save_ifs_2d_fluxes=1

#################################################################################
#
#  PURPOSE: Extract all possible surface flux components of freshwater (FWF) and
#  heat (HTF) from IFS/EC-Earth (native grib output). Calculate their surface
#  integrals over the ocean and land... And output them as monthly time-series
#  into 2 netcdf files.
#
#  Dependencies:
#  CDO, NCO
#
#  Author: L. Broddeau for BaraKuda, November 2016.
#
#
# Surface-integrated Freshwater fluxes over the ocean in Sverdrup (Sv)
# --------------------------------------------------------------------
#
# [ 1Sv==1.E6 m^3/s ]
#
# Lagerloef and colleagues note that:
#
# Evaporation from the global ocean is estimated to be ~ 13 Sv, and the
# precipitation sums to ~ 12.2 Sv. The difference of 0.8 Sv compares with the
# estimate of river input of 1.2 Sv. The apparent imbalance of ~ 0.4 Sv excess
# is smaller than the estimated error bars. Sea level rise due to melting
# glaciers is only ~ 0.01 Sv, so cannot account for the imbalance. Global
# groundwater flows are poorly known but generally estimated to be similarly
# small (Cable et al., 1996). Likely sources of error include scant data in the
# southern oceans and a possible underestimate of evaporation in very-high-wind
# conditions. Other surface flux climatologies display similar patterns, but the
# range of the estimates is quite large and unlikely to be significantly
# improved in the near future.
#
# References: Lagerloef, G., Schmitt, R., Schanze, J. Kao, H-Y, 2010, The Ocean
#             and the Global Water Cycle, Oceanography, Vol.23, No.4
#
#
# Surface integrated Heat fluxes over the ocean in PetaWatts (PW) 
# ---------------------------------------------------------------
#
# [ 1PW==1.E15 Watts ]
# => expect an order of magnitude of 10 PW (Solar ~ 55 PW)
#
##################################################################################

cmsg="ERROR: $0 => global variable"
if [ -z ${EXP} ]; then echo "${cmsg} EXP is unknown!"; exit ; fi
if [ -z ${Y_INI_EC} ]; then echo "${cmsg} Y_INI_EC is unknown!"; exit ; fi
if [ -z ${cyear} ]; then echo "${cmsg} cyear is unknown!"; exit ; fi
if [ -z ${NEMO_OUT_D} ]; then echo "${cmsg} NEMO_OUT_D is unknown!"; exit ; fi
if [ -z ${DIAG_D} ]; then echo "${cmsg} DIAG_D is unknown!"; exit ; fi

echo ; echo "DIAG_D = ${DIAG_D}"; echo
mkdir -p ${DIAG_D}


echo " MOD_CDO => ${MOD_CDO} !!!"
if [ ! "${MOD_CDO}" = "" ]; then module add ${MOD_CDO}; fi

EXP_DIR=`echo ${NEMO_OUT_D} | sed -e "s|/output/nemo||g"`
IFS_OUT_D=`echo ${NEMO_OUT_D} | sed -e "s|/output/nemo|/output/ifs|g"`

FLX_EXTRACT="E,LSP,CP,SSR,STR,SLHF,SSHF"

echo
echo " EXP_DIR = ${EXP_DIR}"
echo " IFS_OUT_D = ${IFS_OUT_D}"
echo " Fields to extract from IFS GG files => ${FLX_EXTRACT}"
echo

YDIR=$((${cyear}-${Y_INI_EC}+1))

dir_ece=`printf "%03d" ${YDIR}`
dir_ece=${EXP_DIR}/output/ifs/${dir_ece}

echo ${dir_ece}

F_AREA=${EXP_DIR}/areas.nc
F_MASK=${EXP_DIR}/masks.nc

echo ${F_AREA} ; ls -l ${F_AREA}
echo ${F_MASK} ; ls -l ${F_MASK}


mkdir -p ./IFS

cd ./IFS/

#rm -f *.grb *.nc *.tmp


NRES_IFS=$(((${TRES_IFS}+1)/2)) ; # ex: (T)255 => (A)128
if [ "${NRES_IFS}" = "80" ]; then NRES_IFS="080"; fi

# Create ifs_area_masked:
#echo
#echo "cdo setmisstoc,0 -ifthen -eqc,0 -selvar,A${NRES_IFS}.msk ${F_MASK} -selvar,A${NRES_IFS}.srf ${F_AREA} metrics.nc"
#cdo setmisstoc,0 -ifthen -eqc,0 -selvar,A${NRES_IFS}.msk ${F_MASK} -selvar,A${NRES_IFS}.srf ${F_AREA} metrics.nc
#echo

ncks -O  -h -v A${NRES_IFS}.msk ${F_MASK} -o metrics.nc
ncrename -h -v A${NRES_IFS}.msk,mask  metrics.nc

ncks -A  -h -v A${NRES_IFS}.srf ${F_AREA} -o metrics.nc
ncrename -h -v A${NRES_IFS}.srf,ifs_area_glob metrics.nc

ncwa -h -O -a y metrics.nc -o metrics.nc # remove y of length !

ncap2 -h -A -s "ifs_area_land=mask*ifs_area_glob"      metrics.nc -o metrics.nc
ncap2 -h -A -s "ifs_area_ocean=(1-mask)*ifs_area_glob" metrics.nc -o metrics.nc

# Checking surface of the ocean and continents to be sure...
ncap2 -h -A -s "srf_glob=ifs_area_glob.total(\$x)*1.E-12" metrics.nc -o metrics.nc
ncatted -h -O -a units,srf_glob,o,c,'10^6 km^2' metrics.nc
ncap2 -h -A -s "srf_ocean=ifs_area_ocean.total(\$x)*1.E-12" metrics.nc -o metrics.nc
ncatted -h -O -a units,srf_ocean,o,c,'10^6 km^2' metrics.nc
ncap2 -h -A -s "srf_land=ifs_area_land.total(\$x)*1.E-12" metrics.nc -o metrics.nc
ncatted -h -O -a units,srf_land,o,c,'10^6 km^2' metrics.nc



# Add degenerate time record:
ncecat   -h -O metrics.nc -o metrics.nc
ncrename -h -d record,time metrics.nc


# CONTROL SECTION FOR IFS:
echo;
srf_glob_ifs=`cdo --no_warnings output -fldsum -selvar,srf_glob metrics.nc 2>/dev/null`
srf_ocean_ifs=`cdo --no_warnings output -fldsum -selvar,srf_ocean metrics.nc 2>/dev/null`
srf_land_ifs=`cdo --no_warnings output -fldsum -selvar,srf_land metrics.nc 2>/dev/null`
echo " *** DEBUG: extract_ifs_surf_fluxes.sh => In IFS data we have:"
echo "            Surface of Earth      = ${srf_glob_ifs} 10^6 km^2" 
echo "            Surface of Ocean      = ${srf_ocean_ifs} 10^6 km^2" 
echo "            Surface of Continents = ${srf_land_ifs} 10^6 km^2" 
echo


wc=`which cdo`
if [ "${wc}" = "" ]; then
    echo "=========================================================="
    echo "ERROR: $0 => No CDO software!!! (command 'cdo' not found)"
    echo "=========================================================="
    echo
    exit
fi

wc=`which ncks`
if [ "${wc}" = "" ]; then
    echo "=========================================================="
    echo "ERROR: $0 => No NCO software!!! (command 'ncks' not found)"
    echo "=========================================================="
    echo
    exit
fi

wc=`which nccopy`
if [ "${wc}" = "" ]; then
    echo "=========================================================="
    echo "WARNING: $0 => command 'nccopy' not found)"
    echo "=========================================================="
    echo
    lnc4=false
else
    lnc4=true
fi



FLX_EXTRACT_S=`echo "${FLX_EXTRACT}" | sed -e s/','/' '/g` ; # with a ' ' instead of a ','



flsm_exp=${EXP_DIR}/output/ifs/001/ICMGG${EXP}+000000
if [ "${i_save_ifs_2d_fluxes}" = "1" ]; then
    
    cdir_ifs_flx=${DIAG_D}/sflx_ifs
    mkdir -p ${cdir_ifs_flx}

    if [ ! -f ${flsm_exp} ]; then
        echo " WARNING: ${0} => you set i_save_ifs_2d_fluxes=1 but file with LSM cannot be found:"
        echo "      ==> ${flsm_exp}"
        echo "      ==> setting i_save_ifs_2d_fluxes=0 !!!"
    else
        echo "cdo -R -t ecmwf -f nc -selvar,LSM ${flsm_exp} LSM_${EXP}_xy.nc"
        cdo -R -t ecmwf -f nc -selvar,LSM ${flsm_exp} LSM_${EXP}_xy.nc
        echo
        # Creating ocean mask MSK:
        ncap2 -h -A -s "MSK=LSM*0."           LSM_${EXP}_xy.nc -o LSM_${EXP}_xy.nc
        ncap2 -O -s 'where(LSM < 0.5) MSK=1.' LSM_${EXP}_xy.nc -o LSM_${EXP}_xy.nc        
    fi
fi



for cm in "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12"; do

    echo
    echo " ${0} => ${cyear}/${cm}"

    fgrb=${dir_ece}/ICMGG${EXP}+${cyear}${cm}

    if [ "${cm}" = "01" ]; then
        pptime=$(cdo showtime -seltimestep,1,2 ${fgrb} | tr -s ' ' ':' | awk -F: '{print ($5-$2)*3600+($6-$3)*60+($7-$4)}' )
        if [ $pptime -le 0 ]; then
            pptime=21600 # default 6-hr output timestep
        fi
        echo " pptime = ${pptime} seconds !"
    fi

    FALL=ALL_${EXP}_${cyear}${cm}

    # Extracting variables of interest and converting to netcdf at the same time (keep gaussian-reduced grid!!!)
    echo "cdo -t ecmwf -f nc -selvar,${FLX_EXTRACT} ${fgrb} ${FALL}.nc"
    cdo -t ecmwf -f nc -selvar,${FLX_EXTRACT} ${fgrb} ${FALL}.nc
    
    if [ "${i_save_ifs_2d_fluxes}" = "1" ]; then
        echo "cdo -R -t ecmwf -f nc -selvar,${FLX_EXTRACT} ${fgrb} ${FALL}_xy.nc"
        cdo -R -t ecmwf -f nc -selvar,${FLX_EXTRACT} ${fgrb} ${FALL}_xy.nc
    fi
    echo

    ncrename -h -O -d rgrid,x ${FALL}.nc

    icpt=0
    for BVAR in ${FLX_EXTRACT_S}; do

        ((icpt++))

        VAR=`echo ${BVAR} | tr '[:upper:]' '[:lower:]'`

        ftreat=${VAR}_${EXP}_${cyear}${cm}

        #BVAR=`echo ${VAR}| tr '[:lower:]' '[:upper:]'`

        # To netcdf monthly:
        #echo "ncra -h -O -v ${BVAR} ${FALL}.nc -O ${ftreat}_m.nc"
        ncra -h -O -v ${BVAR} ${FALL}.nc -o ${ftreat}_m.nc

        # To a flux (originally cumulated flux over time):
        #echo "ncap2 -h -A -s ${VAR}=${BVAR}/${pptime} ${ftreat}_m.nc -o ${ftreat}.nc"
        ncap2 -h -A -s "${VAR}=${BVAR}/${pptime}" ${ftreat}_m.nc -o ${ftreat}.nc ; rm ${ftreat}_m.nc

        # Append ocean mask surface into the file:
        ncks -h -A -v ifs_area_glob  metrics.nc -o ${ftreat}.nc
        ncks -h -A -v ifs_area_ocean metrics.nc -o ${ftreat}.nc
        ncks -h -A -v ifs_area_land  metrics.nc -o ${ftreat}.nc
        

        sign="1."
        if [ "${VAR}" = "e" ]; then sign="-1."; fi

        cu='pw'; cun='PW'; rfact="1.E-15"; ctype="htf" ; cun2d='W/m^2' # ???
        if [ "${VAR}" = "e" ] || [ "${VAR}" = "lsp" ] || [ "${VAR}" = "cp" ]; then
            cu='sv'; cun='Sv'; rfact="1.E-6"; ctype="fwf" ; cun2d='m/s???' # ???
        fi



        # Global horizontal integration
        #------------------------------

        # Multiplying ${VAR} and ifs_area_masked:
        ncap2 -h -A -s "${VAR}2d=(${sign}*ifs_area_ocean*${VAR})"     ${ftreat}.nc -o ${ftreat}.nc
        ncap2 -h -A -s "${VAR}2d_glb=(${sign}*ifs_area_glob*${VAR})"  ${ftreat}.nc -o ${ftreat}.nc
        ncap2 -h -A -s "${VAR}2d_land=(${sign}*ifs_area_land*${VAR})" ${ftreat}.nc -o ${ftreat}.nc

        # Total volume evaporated over ocean during the current month:
        ncap2 -h -A -s "flx_${VAR}_${cu}=${VAR}2d.total(\$x)*${rfact}"           ${ftreat}.nc -o ${ftreat}.nc
        ncap2 -h -A -s "flx_${VAR}_glb_${cu}=${VAR}2d_glb.total(\$x)*${rfact}"   ${ftreat}.nc -o ${ftreat}.nc
        ncap2 -h -A -s "flx_${VAR}_land_${cu}=${VAR}2d_land.total(\$x)*${rfact}" ${ftreat}.nc -o ${ftreat}.nc
        ncatted -O -a units,flx_${VAR}_${cu},o,c,"${cun}"     ${ftreat}.nc
        ncatted -O -a units,flx_${VAR}_glb_${cu},o,c,"${cun}" ${ftreat}.nc
        ncatted -O -a units,flx_${VAR}_land_${cu},o,c,"${cun}" ${ftreat}.nc

        ncks -h -A -v flx_${VAR}_${cu}      ${ftreat}.nc -o final_${ctype}_${cm}.nc
        ncks -h -A -v flx_${VAR}_glb_${cu}  ${ftreat}.nc -o final_${ctype}_${cm}.nc
        ncks -h -A -v flx_${VAR}_land_${cu} ${ftreat}.nc -o final_${ctype}_${cm}.nc



        # Same for field in lat-lon:
        if [ "${i_save_ifs_2d_fluxes}" = "1" ]; then
            ncra  -h -O -v ${BVAR} ${FALL}_xy.nc -o ${ftreat}_m_xy.nc ; # Montly !
            ncap2 -h -A -s "tmp=${BVAR}/${pptime}" ${ftreat}_m_xy.nc -o ${ftreat}_xy.nc ; rm ${ftreat}_m_xy.nc
            ncks  -h -A -v MSK LSM_${EXP}_xy.nc  -o ${ftreat}_xy.nc ; # Appends 2D ocean mask !
            # Masking and removing the rif-raf:
            ncap2 -h -A -s "tmpd=${sign}*tmp*MSK-(1.-MSK)*9999." ${ftreat}_xy.nc -o tmp.nc ; rm -f ${ftreat}_xy.nc
            ncap2 -h -A -s "${VAR}=float(tmpd)" tmp.nc -o tmp.nc ; # to float
            ncks  -h -O -v ${VAR} tmp.nc -o ${ftreat}_xy.nc ; rm -f tmp.nc
            ncatted -O -h -a _FillValue,${VAR},o,float,-9999. ${ftreat}_xy.nc
            ncatted -O -a units,${VAR},o,c,"${cun2d}" ${ftreat}_xy.nc ; # correct unit!
        fi

        rm -f ${ftreat}.nc

        # End loop variables
    done

    rm -f ${FALL}.nc

    echo " *** month ${cm} done!"
    echo

done


if [ "${i_save_ifs_2d_fluxes}" = "1" ]; then
    fo2d=FLX_IFS_2D_${EXP}_${cyear}.nc ;    rm -f ${fo2d}
    for BVAR in ${FLX_EXTRACT_S}; do
        VAR=`echo ${BVAR} | tr '[:upper:]' '[:lower:]'`
        echo "ncrcat -h -A ${VAR}_${EXP}_${cyear}*_xy.nc -o ${fo2d}"
        ncrcat -h -A ${VAR}_${EXP}_${cyear}*_xy.nc -o ${fo2d}
        echo
    done
    
    # Cleaning attributes:
    for ga in "history" "NCO" "CDO" "nco_openmp_thread_number" "history_of_appended_files" "CDI" "institution"; do
        ncatted -h -O -a ${ga},global,d,c, ${fo2d}
    done
    ncatted -h -O -a About,global,o,c,"Created by Barakuda (https://github.com/brodeau/barakuda)" ${fo2d}

    if ${lnc4}; then
        nccopy -k 4 -d 9 ${fo2d} ${fo2d}4 &
    fi
    rm -f *_${EXP}_${cyear}*_xy.nc
fi



for ct in "htf" "fwf"; do
    echo "ncrcat -O final_${ct}_*.nc -o final_${ct}.nc"
    ncrcat -O final_${ct}_*.nc -o final_${ct}.nc

    ncap2 -h -O -s "time=array(${cyear}.0416667,0.08333333,\$time)" final_${ct}.nc -o final_${ct}.nc
    ncatted -O -a units,time,o,c,'years' final_${ct}.nc
done


# Freshwater fluxes that need to be built:
ncap2 -h -A -s "flx_p_sv=flx_cp_sv+flx_lsp_sv" final_fwf.nc  ; # total precip over ocean
ncap2 -h -A -s "flx_emp_sv=flx_e_sv-flx_p_sv"  final_fwf.nc  ; # E-P over ocean
# Same, Ocean+Land:
ncap2 -h -A -s "flx_p_glb_sv=flx_cp_glb_sv+flx_lsp_glb_sv" final_fwf.nc
ncap2 -h -A -s "flx_emp_glb_sv=flx_e_glb_sv-flx_p_glb_sv"  final_fwf.nc
# Same, Land:
ncap2 -h -A -s "flx_p_land_sv=flx_cp_land_sv+flx_lsp_land_sv" final_fwf.nc
ncap2 -h -A -s "flx_emp_land_sv=flx_e_land_sv-flx_p_land_sv"  final_fwf.nc


# Heat fluxes that need to be built:
ncap2 -h -A -s "flx_qnet_pw=flx_ssr_pw+flx_str_pw+flx_slhf_pw+flx_sshf_pw" final_htf.nc ; # Net heat flux over ocean
# Same, Ocean+Land:
ncap2 -h -A -s "flx_qnet_glb_pw=flx_ssr_glb_pw+flx_str_glb_pw+flx_slhf_glb_pw+flx_sshf_glb_pw" final_htf.nc
# Same, Land:
ncap2 -h -A -s "flx_qnet_land_pw=flx_ssr_land_pw+flx_str_land_pw+flx_slhf_land_pw+flx_sshf_land_pw" final_htf.nc ; # Net heat flux over ocean

rm -f metrics.nc final_htf_*.nc final_fwf_*.nc



for ct in "htf" "fwf"; do

    fout=${DIAG_D}/mean_${ct}_IFS_${EXP}_GLO.nc
    
    # Cleaning attributes:
    for ga in "history" "NCO" "CDO" "nco_openmp_thread_number" "history_of_appended_files"; do
        ncatted -h -O -a ${ga},global,d,c, final_${ct}.nc
    done
    
    if [ ! -f ${fout} ]; then
        ncatted -h -O -a Global_Ocean_Area,global,o,c,"${srf_ocean_ifs} 10^6 km^2" final_${ct}.nc
        ncatted -h -O -a Global_Land_Area,global,o,c,"${srf_land_ifs} 10^6 km^2"   final_${ct}.nc
        ncatted -h -O -a Global_Area,global,o,c,"${srf_glob_ifs} 10^6 km^2"        final_${ct}.nc
        ncatted -h -O -a About,global,o,c,"Created by Barakuda (https://github.com/brodeau/barakuda)" final_${ct}.nc
        mv final_${ct}.nc ${fout}
    else
        ncrcat -h -A ${fout} final_${ct}.nc -o ${fout}
    fi

    rm -f final_${ct}.nc

done


if [ "${i_save_ifs_2d_fluxes}" = "1" ]; then
    if ${lnc4}; then
        wait
        rm -f ${fo2d}
    fi
    mv -f ${fo2d}* ${cdir_ifs_flx}/
fi

if [ ! "${MOD_CDO}" = "" ]; then module rm ${MOD_CDO}; fi

cd ../
rm -rf ./IFS

exit 0
