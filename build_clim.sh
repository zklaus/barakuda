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

ivt=1   ; # Create a climatology for VT
iamoc=1 ; # Create a climatology for 2D lat-depth AMOC?
ibpsi=0 ; # Create a climatology for barotropic stream function

export script=build_clim
#export BARAKUDA_ROOT=`pwd`
[ -z ${BARAKUDA_ROOT+x} ] && export BARAKUDA_ROOT=${PWD}

# Checking available configs
list_conf=`\ls configs/config_*.sh` ; list_conf=`echo ${list_conf} | sed -e s/'configs\/config_'/''/g -e s/'.sh'/''/g`

# Important bash functions:
. ${BARAKUDA_ROOT}/src/bash/bash_functions.bash

barakuda_init

while getopts R:f:i:e:C:h option ; do
    case $option in
        R) export EXP=${OPTARG};;
        f) export IFREQ_SAV_YEARS=${OPTARG} ;;
        i) export Y1=${OPTARG} ;;
        e) export Y2=${OPTARG} ;;
        C) export CONFIG=${OPTARG};;
        h)  usage;;
        \?) usage ;;
    esac
done

barakuda_check

# sourcing configuration file
fconfig=${BARAKUDA_ROOT}/configs/config_${CONFIG}.sh

if [ -f ${fconfig} ]; then
    echo "Sourcing ${fconfig} !"
    . ${fconfig}
else
    echo "PROBLEM: cannot find file ${fconfig} !"; exit
fi

# If auto-submit experiment (ece_exp=10) then overides a few functions with:
if [ ${ece_exp} -ge 10 ]; then
    . ${BARAKUDA_ROOT}/src/bash/bash_functions_autosub.bash
fi

# If 3D fieds are annual averaged then overides a few functions with:
if [ ! "${ANNUAL_3D}" = "" ]; then
    . ${BARAKUDA_ROOT}/src/bash/bash_functions_1y.bash
fi

barakuda_setup

CPRMN="${CPREF}"
CPRAN=""
if [ ! "${ANNUAL_3D}" = "" ]; then
    CPRAN=`echo ${CPREF} | sed -e s/"_${TSTAMP}_"/"_${ANNUAL_3D}_"/g`
fi

echo
echo " SETTINGS: "
echo "   *** DIAG_D   = ${DIAG_D} "
echo "   *** CLIM_DIR = ${CLIM_DIR} "
echo "   *** TMP_DIR  = ${TMP_DIR} "
echo "   *** Y1       = ${Y1} "
echo "   *** Y2       = ${Y2} "
echo "   *** CONFIG   = ${CONFIG} "
echo "   *** GRID     = ${ORCA} "
echo "   *** CONFEXP  = ${CONFEXP} "
echo "   *** EXP      = ${EXP} "
echo "   *** CPREF    = ${CPREF} "
echo "   *** Monthly files prefix = ${CPRMN} (${NEMO_SAVED_FILES})"
if [ ! "${CPRAN}" = "" ]; then
    echo "   ***  Annual files prefix = ${CPRAN} (${NEMO_SAVED_FILES_3D})"
fi
echo "   *** IFREQ_SAV_YEARS = ${IFREQ_SAV_YEARS} "
echo "   *** NCDF_DIR        = ${NCDF_DIR} "
echo

wc=`which nccopy`
if [ "${wc}" = "" ]; then
    CP2NC4="${NCDF_DIR}/bin/nccopy -k 4 -d 9"
else
    CP2NC4="nccopy -k 4 -d 9"
fi

# cdftools execs are found:
export PATH=${BARAKUDA_ROOT}/cdftools_light/bin:${PATH}


Y1=$((${Y1}+0))
Y2=$((${Y2}+0))
CY1=`printf "%04d" ${Y1}`
CY2=`printf "%04d" ${Y2}`


mkdir -p ${CLIM_DIR}


# Variables to extract:
V2E="${NN_SST},${NN_SSS},${NN_SSH},${NN_T},${NN_S},${NN_MLD}"
C2ET="nav_lon,nav_lat,deptht"
C2EU="nav_lon,nav_lat,depthu"
C2EV="nav_lon,nav_lat,depthv"
C2EW="nav_lon,nav_lat,depthw"

GRID_IMP="grid_T"
if [ ${ivt} -eq 1 ] || [ ${ibpsi} -eq 1 ]; then
    GRID_IMP+=" grid_U"
fi
if [ ${iamoc} -eq 1 ] || [ ${ivt} -eq 1 ] || [ ${ibpsi} -eq 1 ]; then
    GRID_IMP+=" grid_V"
fi
if [ `contains_string ${FILE_ICE_SUFFIX} ${NEMO_SAVED_FILES}` -eq 1 ]; then
    GRID_IMP+=" ${FILE_ICE_SUFFIX}"
fi
if [ `contains_string SBC ${NEMO_SAVED_FILES}` -eq 1 ]; then
    GRID_IMP+=" SBC"
fi
echo; echo " GRID_IMP = ${GRID_IMP}"; echo


# Checking what files we have / plan to use:
if [ -z "${NEMO_SAVED_FILES}" ]; then
    echo "Please specify which NEMO files are saved (file suffixes, grid_T, ..., icemod) ?"
    echo " => set the variable NEMO_SAVED_FILES in your config_${CONFIG}.sh file!"; exit
fi
VAF=( "grid_T" "grid_U" "grid_V" "icemod" "SBC" )
js=0 ; gimp_new=""
for sf in ${VAF[*]}; do
    echo "Checking ${sf}..."
    ca=`echo "${NEMO_SAVED_FILES} ${NEMO_SAVED_FILES_3D}" | grep ${sf}`
    cb=`echo "${GRID_IMP}"         | grep ${sf}`
    if [ "${ca}" = "" ]; then
        if [ "${cb}" != "" ]; then
            echo "PROBLEM! The diags you specified say you need ${sf} files"
            echo "     => but you have not specified ${sf} in NEMO_SAVED_FILES !"; exit
        fi
    else
        gimp_new="${sf} ${gimp_new}"
    fi
    ((js++))
done
GRID_IMP=${gimp_new}
echo; echo "File types to import: ${GRID_IMP}"; echo; echo


VCM=( "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" )

cd ${TMP_DIR}/

echo; echo "In:"; pwd

barakuda_import_mesh_mask

ls -l ; echo; echo


nby=$((${Y2}-${Y1}+1))

if [ ! $((${nby}%${IFREQ_SAV_YEARS})) -eq 0 ]; then
    echo " Number of years should be a multiple of ${IFREQ_SAV_YEARS}!"; exit
fi

if [ ${ece_exp} -gt 0 ]; then
    if [ ! -d ${NEMO_OUT_D}/001 ]; then echo "ERROR: since ece_exp=${ece_exp}, there should be a directory 001 in:"; echo " ${NEMO_OUT_D}"; fi
fi

ffirsty="${DIAG_D}/first_year.info"
if [ ! -f  ${ffirsty} ]; then echo "ERROR: file ${ffirsty} not found!!!"; exit; fi
export YEAR_INI_F=`cat ${ffirsty}`
export jyear=${Y1}


while [ ${jyear} -le ${Y2} ]; do

    export cyear=`printf "%04d" ${jyear}` ; echo ; echo "Year = ${cyear}"

    cpf=""
    if [ ${ece_exp} -gt 0 ]; then
        iy=$((${jyear}-${Y1}+1+${Y1}-${YEAR_INI_F}))
        dir_ece=`printf "%03d" ${iy}`
        echo " *** ${cyear} => dir_ece = ${dir_ece}"
        cpf="${dir_ece}/"
    fi

    TTAG_ann=${cyear}0101_${cyear}1231

    i_get_file=0
    if [ $((${jyear}%${IFREQ_SAV_YEARS})) -eq 0 ]; then
        barakuda_check_year_is_complete  ; # lcontinue might be updated to false!
    fi

    CRTM=${CPRMN}${TTAG}
    CRT1M=${CPRMN}${TTAG_ann}

    barakuda_import_files

    # Monthly files to work with for current year:
    ft1m=${CRT1M}_grid_T.nc
    fu1m=${CRT1M}_grid_U.nc
    fv1m=${CRT1M}_grid_V.nc    
    # Annual files to work with for current year:
    CRT1Y=`echo ${CRT1M} | sed -e s/"_${TSTAMP}_"/"_${ANNUAL_3D}_"/g`
    ft1y=${CRT1Y}_grid_T.nc
    fu1y=${CRT1Y}_grid_U.nc
    fv1y=${CRT1Y}_grid_V.nc
    fj1y=${CRT1Y}_${FILE_ICE_SUFFIX}.nc ; # can be icemod or grid_T ....
    CFG3D=${CRT1M}
    CPREF3D=${CPRMN}
    #
    # Files that contain the 3D fields (might be monthly or annaual sometimes (when "${ANNUAL_3D}" = "1y")    $
    ft3d=${ft1m}
    fu3d=${fu1m}
    fv3d=${fv1m}
    if [ "${ANNUAL_3D}" = "1y" ]; then
        [[ ${NEMO_SAVED_FILES_3D} =~ (^|[[:space:]])"grid_U"($|[[:space:]]) ]] \
            && CPREF3D=${CPRAN}; CFG3D=${CRT1Y}; ft3d=${ft1y}; fu3d=${fu1y}; fv3d=${fv1y} \
            || echo "...default"
        echo ""
    fi
    fvt=${CFG3D}_VT.nc

    echo

    if [ ${ivt} -eq 1 ]; then
        # Creating VT files:
        echo " *** CALLING: cdfvT.x ${CFG3D} ${NN_T} ${NN_S} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV} &"
        cdfvT.x ${CFG3D} ${NN_T} ${NN_S} ${NN_U} ${NN_V} ${NN_U_EIV} ${NN_V_EIV} &
        echo
    fi

    echo

    if [ ${iamoc} -eq 1 ]; then
        rm -f moc.nc
        echo " *** CALLING: cdfmoc.x ${fv3d} ${NN_V} ${NN_V_EIV} &"
        cdfmoc.x ${fv3d} ${NN_V} ${NN_V_EIV} &
        echo
    fi

    if [ ${ibpsi} -eq 1 ]; then
        rm -f psi.nc
        echo " *** CALLING: cdfpsi.x ${fu3d} ${fv3d} ${NN_U} ${NN_V} V &"
        cdfpsi.x ${fu3d} ${fv3d} ${NN_U} ${NN_V} V &
        echo
    fi

    wait

    ncwa -O -a x moc.nc ${CFG3D}_MOC.nc ; # removing degenerate x record...
    rm -f moc.nc
    
    if [ ${ibpsi} -eq 1 ]; then mv -f psi.nc ${CFG3D}_PSI.nc; fi
    
    echo "After year ${jyear}:"; ls -l *.nc*
    echo
    
    ((jyear++))
    export jyear=${jyear}

done

echo
echo "Phase 2:"; ls ; echo


# Mean monthly climatology

SUFF_FOR_MONTHLY="${NEMO_SAVED_FILES}"
if [ "${ANNUAL_3D}" = "" ]; then SUFF_FOR_MONTHLY+=" VT MOC PSI"; fi
echo; echo "Will build monthly clim for files with these suffixes:"; echo ${SUFF_FOR_MONTHLY}; echo

for suff in ${SUFF_FOR_MONTHLY}; do
    if [ -f ./${CRT1M}_${suff}.nc ]; then
        echo ; echo " Treating ${suff} files!"; echo
        f2c=mclim_${CONFEXP}_${CY1}-${CY2}_${suff}.nc
        rm -f ${CLIM_DIR}/${f2c}*
        echo
        jm=0
        for cm in ${VCM[*]}; do
            ((jm++))
            if [ -f ./${CRT1M}_${suff}.nc ]; then
                echo; ls ; echo
                echo "ncra -F -O -d time_counter,${jm},,12 ${CPRMN}*0101_*1231_${suff}.nc -o mean_m${cm}_${suff}.nc"
                ncra -F -O -d time_counter,${jm},,12 ${CPRMN}*0101_*1231_${suff}.nc -o mean_m${cm}_${suff}.nc &
                echo
            fi
        done
        wait

        echo; ls ; echo
        echo "ncrcat -O  mean_m*_${suff}.nc out_${suff}.nc"
        ncrcat -O  mean_m*_${suff}.nc out_${suff}.nc
        rm mean_m*_${suff}.nc
        echo
        echo "mv -f out_${suff}.nc ${f2c}"
        mv -f out_${suff}.nc ${f2c}
        echo
    else
        echo ; echo " Ignoring monthly ${suff} files!"; echo
    fi
    
done ; # loop along monthly files suffixes

wait


# Mean annual climatology for 3D-variable-based fields:

SUFF_FOR_ANNUAL="VT MOC PSI"
if [ ! "${ANNUAL_3D}" = "" ]; then SUFF_FOR_ANNUAL+=" ${NEMO_SAVED_FILES_3D}"; fi
echo; echo "Will build annual clim for files with these suffixes:"; echo ${SUFF_FOR_ANNUAL}; echo

for suff in ${SUFF_FOR_ANNUAL}; do
    if [ -f ./${CRT1Y}_${suff}.nc ] || [ -f ./${CRT1M}_${suff}.nc ]; then
        echo ; echo " Treating ${suff} files!"; echo
        f2c=aclim_${CONFEXP}_${CY1}-${CY2}_${suff}.nc
        rm -f ${CLIM_DIR}/${f2c}*
        #
        echo "ncra -O ${CPREF3D}*_${suff}.nc -o ${f2c} &"
        ncra -O ${CPREF3D}*_${suff}.nc -o ${f2c} &
        echo
    else
        echo ; echo " Ignoring annual ${suff} files!"; echo
    fi    
done ; # loop along annual files suffixes
wait

rm -f ${CPRMN}*0101_*1231_*.nc ${CPRAN}*0101_*1231_*.nc



list=`\ls [am]clim_${CONFEXP}*.nc`
for ff in ${list}; do
    fn=`echo ${ff} | sed -e s/'.nc'/'.nc4'/g`
    echo "${CP2NC4}  ${ff} ${fn} &"
    ${CP2NC4}  ${ff} ${fn} &
    echo
done

wait

for cl in aclim mclim; do
    echo "mv -f ${cl}_${CONFEXP}*.nc4 ${CLIM_DIR}/"
    mv -f ${cl}_${CONFEXP}*.nc4 ${CLIM_DIR}/
    echo
done


echo;echo
echo "${CY1}-${CY2}" > ${CLIM_DIR}/last_clim
echo "Climatology saved into: ${CLIM_DIR}/"
echo;echo

cd /tmp/
rm -rf ${TMP_DIR} 2>/dev/null

exit 0
