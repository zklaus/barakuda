#!/bin/bash

#  B E T A  ! ! !

# Diag to test:
# One at a time please!!!
itrsig=0
imhst=0
iamoc=0
icrosssect=0
itempsal=0
ifsflx=0
imean2d=1
imean3d=0
iSflx=0
ienso=0
imov=0
issh=0
iwind=0
its=0
imld=0
irnf=0
iice=0
iemp=0
icmip5=0
ihov=0

#CONFIG="ORCA1_L75"
#ARCH="T159_ece32_marenostrum"
#export EXP="LR1E" ; NC=nc4 ; jyear=1990

#CONFIG="ORCA1_L42"
#ARCH="ece22_triolith"
#export EXP="SPIN" ; NC=nc4 ; jyear=2540

#CONFIG="ORCA2_L31"
#ARCH="ece32_marenostrum"
#export EXP="LR20" ; NC=nc4 ; jyear=2540

#CONFIG="ORCA025_L75"
#ARCH="T511_ece32_triolith"
#export EXP="HC71" ; NC=nc ; jyear=1990

#CONFIG="ORCA025_L75"
#ARCH="etienne"
#export EXP="a0ez" ; NC=nc ; jyear=1945

#CONFIG="ORCA1_L75"
#ARCH="ro10"
#export EXP="ro10" ; NC=nc ; jyear=1995

CONFIG="ORCA1_L75"
ARCH="uwe"
export EXP="SN71" ; NC=nc ; jyear=1995

#CONFIG="ORCA1_L75"
#ARCH="T159_ece32_triolith"
#export EXP="LB30" ; NC=nc4 ; jyear=2010



export BARAKUDA_ROOT=`pwd | sed -e "s|/python||g"`

. ${BARAKUDA_ROOT}/src/bash/bash_functions.bash
. ${BARAKUDA_ROOT}/configs/config_${CONFIG}_${ARCH}.sh

ORCA_LIST="ORCA025.L75 ORCA1.L75 ORCA1.L46 ORCA1.L42 ORCA2.L31 ORCA2.L46"

for og in ${ORCA_LIST}; do
    ca=""; ca=`echo ${CONFIG} | grep ${og}` ; if [ "${ca}" != "" ]; then export ORCA=${og}; fi
done
if [ "${ORCA}" = "" ]; then echo "ORCA grid of config ${CONFIG} not supported yet"; exit; fi

export CONFEXP=${ORCA}-${EXP}
export DIAG_D=${DIAG_DIR}/${CONFEXP} ; mkdir -p ${DIAG_D}

echo ; echo " *** DIAG_D = ${DIAG_D} !"; echo

HERE=`pwd`

finfoclim=${DIAG_D}/clim/last_clim

y1_clim=`cat ${finfoclim} | cut -d - -f1`
y2_clim=`cat ${finfoclim} | cut -d - -f2`

export COMP2D="OBS"

# To know the name of NEMO output files:
export NEMO_OUT_D=`echo ${NEMO_OUT_STRCT} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<EXP>|${EXP}|g"`
if [ ! -d ${NEMO_OUT_D} ]; then echo "Unfortunately we could not find ${NEMO_OUT_D}"; exit; fi
YEAR_INI=1990 ; YEAR_INI_F=1990
export cyear=`printf "%04d" ${jyear}`
if [ ${ece_exp} -gt 0 ]; then
    iy=$((${jyear}-${YEAR_INI}+1+${YEAR_INI}-${YEAR_INI_F}))
    dir_ece="`printf "%03d" ${iy}`/"
fi
CPREF=`echo ${NEMO_FILE_PREFIX} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<EXP>|${EXP}|g" -e "s|<TSTAMP>|${TSTAMP}|g"`

if [ ${icrosssect} -eq 1 ] || [ ${imean2d} -eq 1 ] || [ ${imov} -eq 1 ]; then
    ft=${NEMO_OUT_D}/${dir_ece}${CPREF}${cyear}0101_${cyear}1231_grid_T.${NC}
    check_if_file ${ft}
    fj=${NEMO_OUT_D}/${dir_ece}${CPREF}${cyear}0101_${cyear}1231_${FILE_ICE_SUFFIX}.${NC}
    check_if_file ${fj}
fi

if [ ${imean3d} -eq 1 ]; then
    ft=${NEMO_OUT_D}/${dir_ece}${CPREF}${cyear}0101_${cyear}1231_grid_T.${NC}
    check_if_file ${ft}
fi

if [ ${iSflx} -eq 1 ]; then
    ft=${NEMO_OUT_D}/${dir_ece}${CPREF}${cyear}0101_${cyear}1231_${FILE_FLX_SUFFIX}.${NC}
    check_if_file ${ft}
fi



export PYTHONPATH=${PYTHON_HOME}/lib/python2.7/site-packages:${BARAKUDA_ROOT}/python/modules

echo ; echo " *** DIAG_D=${DIAG_D} !"; echo


rm -f *.png

# Time for diags:

if [ ${ienso} -eq 1 ]; then
    CMD="python exec/plot_enso.py ${DIAG_D}/Nino34_${CONFEXP}.nc ${NN_SST}"
    echo ; echo " CMD = ${CMD} "; echo
fi

if [ ${imhst} -eq 1 ]; then
    CMD="python exec/plot_hovm_merid_trsp.py"
    echo ; echo " CMD = ${CMD} "; echo
fi

if [ ${iamoc} -eq 1 ]; then
    if [ ! -f ${DIAG_D}/clim/last_clim ]; then echo "Boooo!"; exit; fi
    CLIM_PER=`cat ${DIAG_D}/clim/last_clim`
    iclyear=`echo ${CLIM_PER} | sed -e s/'-'/' '/g`
    CMD="python exec/moc.py ${iclyear}"
    echo ; echo " CMD = ${CMD} "; echo
fi

if [ ${icrosssect} -eq 1 ]; then
    export DIAG_D=`pwd`
    CMD="python exec/cross_sections.py ${ft} ${jyear}"
fi


if [ ${itempsal} -eq 1 ]; then
    if [ ! -f ${DIAG_D}/clim/last_clim ]; then echo "Boooo!"; exit; fi
    CLIM_PER=`cat ${DIAG_D}/clim/last_clim`
    iclyear=`echo ${CLIM_PER} | sed -e s/'-'/' '/g`
    CMD="python exec/temp_sal.py ${iclyear}"
    echo ; echo " CMD = ${CMD} "; echo
fi

if [ ${its} -eq 1 ]; then
    #diag=3d_thetao ; ln -sf ${DIAG_D}/3d_${NN_T}*.nc .
    #diag=mean_zos  ; ln -sf ${DIAG_D}/mean_${NN_SSH}*.nc .
    #diag=mean_htf ; ln -sf ${DIAG_D}/mean_htf*.nc .
    #diag=mean_fwf ; ln -sf ${DIAG_D}/mean_fwf*.nc .
    diag=transport_sections ; ln -sf ${DIAG_D}/transport_*sect_*.nc .
    #ln -sf ${DIAG_D}/${diag}*.nc .
    CMD="python exec/plot_time_series.py ${diag}"
fi

if [ ${itrsig} -eq 1 ]; then
    CMD="python exec/plot_trsp_sigma.py"
fi



if [ ${ifsflx} -eq 1 ]; then
    export DIAG_D="`pwd`/flx"
    CMD="${BARAKUDA_ROOT}/src/bash/extract_ifs_surf_fluxes.sh"
fi

if [ ${imean2d} -eq 1 ]; then
    export DIAG_D=`pwd`
    CMD="python exec/mean_2d.py ${ft} ${jyear}"
fi

if [ ${imean3d} -eq 1 ]; then
    export DIAG_D=`pwd`
    CMD="python exec/mean_3d.py ${ft} ${jyear} T"
fi

if [ ${iSflx} -eq 1 ]; then
    export DIAG_D=`pwd`
    CMD="python exec/flux_int_basins.py ${ft} ${jyear}"
fi


if [ ${imov} -eq 1 ]; then
    #for cv in sst mld sss; do
    cv="ice"
    #python exec/prepare_movies.py ${ft} ${jyear} ${cv}
    python ./prepare_movies.py ${ft} ${jyear} ${cv}
    #done
    exit
fi

if [ ${issh} -eq 1 ]; then
    CMD="python exec/ssh.py ${y1_clim} ${y2_clim}"
fi

if [ ${iwind} -eq 1 ]; then
    CMD="python exec/wind.py ${y1_clim} ${y2_clim}"
fi


if [ ${imld} -eq 1 ]; then
    CMD="python exec/mld.py ${y1_clim} ${y2_clim}"
fi



echo
echo "DOING: ${CMD}"
${CMD}


# Add other diags here:






exit
# BELOW = OLD STUFFS, fix!










if [ ${ihov} -eq 1 ]; then
    export EXP=cp70
    export ORCA=ORCA1.L75
    export DIAG_D=/proj/bolinc/users/x_laubr/tmp/barakuda/ORCA1.L75_ece32b/ORCA1.L75-${EXP}
    export MM_FILE=/proj/bolinc/users/x_laubr/klaus/mesh_mask.nc
    export BM_FILE=/proj/bolinc/users/x_laubr/klaus/basin_mask.nc
    export NN_T="thetao"
    export NN_S="so"
    #
    cd ${DIAG_D}/
    python /home/x_laubr/DEV/barakuda/python/exec/plot_hovm_tz.py 1996 2000

    mv -f hov_*_ORCA1.L75-${EXP}*.png ${HERE}/
    #
fi








if [ ${iemp} -eq 1 ]; then

    export ORCA="ORCA1.L75"
    #export EXP="32bI"
    export EXP="cp00"
    export TSTAMP="1m"
    export DIAG_D="/proj/bolinc/users/x_laubr/tmp/barakuda/ORCA1.L75_ece32b/ORCA1.L75-${EXP}"
    #export NN_RNF="runoffs"
    export MM_FILE="/proj/bolinc/users/x_laubr/tmp/barakuda/test/mesh_mask.nc"
    export TRANSPORT_SECTION_FILE="boo"
    export LMOCLAT="boo" ; export NN_SSH="boo" ; export NN_SSS="boo" ; export NN_S="boo"
    export NN_MLD="boo" ; export NN_SST="boo" ; export NN_T="boo"
    export NN_FWF="wfo"       ; # name of net freshwater flux (E-P-R) in "FILE_FLX_SUFFIX" file...
    export NN_EMP="emp_oce"   ; # name of E-P in "FILE_FLX_SUFFIX" file...
    export NN_P="precip"   ; # name of P in "FILE_FLX_SUFFIX" file...
    export NN_RNF="XXX"   ; # name of continental runoffs in "FILE_FLX_SUFFIX" file...



    export FILE_DEF_BOXES="/home/x_laubr/DEV/barakuda/data/def_boxes_convection_ORCA1.txt"

    cd ${DIAG_D}/

    python /home/x_laubr/DEV/barakuda/python/exec/plot_time_series.py mean_fwf





fi











if [ ${irnf} -eq 1 ]; then

    export ORCA="ORCA1.L75"
    export EXP="LB03"
    export DIAG_D="/proj/bolinc/users/x_laubr/tmp/barakuda/ORCA1.L75_ece32b/ORCA1.L75-${EXP}"
    export NN_RNF="runoffs"
    export MM_FILE="/proj/bolinc/users/x_laubr/tmp/barakuda/test/mesh_mask.nc"

    python exec/runoffs.py 1997 1999

fi




