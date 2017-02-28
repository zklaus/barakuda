#!/bin/bash

#  B E T A  ! ! !

# Diag to test:
imhst=0
iamoc=0
icrosssect=0
itempsal=0
ifsflx=0
imean2d=0
imean3d=0
ienso=0
imov=0
issh=0
its=1
imld=0
irnf=0
iice=0
iemp=0
icmip5=0
ihov=0

#CONFIG="ORCA1_L75"
#ARCH="T159_ece32_marenostrum"
#export RUN="LR1E" ; NC=nc4 ; jyear=1990

#CONFIG="ORCA1_L42"
#ARCH="ece22_triolith"
#export RUN="SPIN" ; NC=nc4 ; jyear=2540

#CONFIG="ORCA2_L31"
#ARCH="ece32_marenostrum"
#export RUN="LR20" ; NC=nc4 ; jyear=2540


#CONFIG="ORCA1_L75"
#ARCH="T159_ece32_triolith"
#export RUN="LB2E" ; NC=nc4 ; jyear=2010

#CONFIG="ORCA025_L75"
#ARCH="T511_ece32_triolith"
#export RUN="HC71" ; NC=nc ; jyear=1990

#CONFIG="ORCA025_L75"
#ARCH="etienne"
#export RUN="a0ez" ; NC=nc ; jyear=1945

CONFIG="ORCA1_L75"
ARCH="T255_ece32_triolith"
export RUN="LB10" ; NC=nc4 ; jyear=1990


export BARAKUDA_ROOT=`pwd | sed -e "s|/python||g"`

. ${BARAKUDA_ROOT}/src/bash/bash_functions.bash
. ${BARAKUDA_ROOT}/configs/config_${CONFIG}_${ARCH}.sh

ORCA_LIST="ORCA025.L75 ORCA1.L75 ORCA1.L46 ORCA1.L42 ORCA2.L31 ORCA2.L46"

for og in ${ORCA_LIST}; do
    ca=""; ca=`echo ${CONFIG} | grep ${og}` ; if [ "${ca}" != "" ]; then export ORCA=${og}; fi
done
if [ "${ORCA}" = "" ]; then echo "ORCA grid of config ${CONFIG} not supported yet"; exit; fi

export CONFRUN=${ORCA}-${RUN}
export DIAG_D=${DIAG_DIR}/${CONFRUN} ; mkdir -p ${DIAG_D}

echo ; echo " *** DIAG_D = ${DIAG_D} !"; echo

HERE=`pwd`

finfoclim=${DIAG_D}/clim/last_clim

y1_clim=`cat ${finfoclim} | cut -d - -f1`
y2_clim=`cat ${finfoclim} | cut -d - -f2`

export COMP2D="CLIM"

# To know the name of NEMO output files:
export NEMO_OUT_D=`echo ${NEMO_OUT_STRCT} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<RUN>|${RUN}|g"`
if [ ! -d ${NEMO_OUT_D} ]; then echo "Unfortunately we could not find ${NEMO_OUT_D}"; exit; fi
YEAR_INI=1990 ; YEAR_INI_F=1990
export cyear=`printf "%04d" ${jyear}`
if [ ${ece_run} -gt 0 ]; then
    iy=$((${jyear}-${YEAR_INI}+1+${YEAR_INI}-${YEAR_INI_F}))
    dir_ece="`printf "%03d" ${iy}`/"
fi
CPREF=`echo ${NEMO_FILE_PREFIX} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<RUN>|${RUN}|g" -e "s|<TSTAMP>|${TSTAMP}|g"`

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



export PYTHONPATH=${PYTHON_HOME}/lib/python2.7/site-packages:${BARAKUDA_ROOT}/python/modules

echo ; echo " *** DIAG_D=${DIAG_D} !"; echo


rm -f *.png

# Time for diags:

if [ ${ienso} -eq 1 ]; then
    CMD="python exec/plot_enso.py ${DIAG_D}/Nino34_${CONFRUN}.nc ${NN_SST}"
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
    export DIAG_D="."
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
    #diag=3d_thetao
    diag=mean_fwf
    #diag=mean_htf
    ln -sf ${DIAG_D}/${diag}*.nc .
    CMD="python exec/plot_time_series.py ${diag}"
fi


if [ ${ifsflx} -eq 1 ]; then
    export DIAG_D="."
    CMD="${BARAKUDA_ROOT}/src/bash/extract_ifs_surf_fluxes.sh"
fi

if [ ${imean2d} -eq 1 ]; then
    export DIAG_D="."
    CMD="python exec/mean_2d.py ${ft} ${jyear}"
fi

if [ ${imean3d} -eq 1 ]; then
    export DIAG_D="."
    CMD="python exec/mean_3d.py ${ft} ${jyear} T"
fi


if [ ${imov} -eq 1 ]; then
    for cv in ice sst sss; do
        python exec/prepare_movies.py ${ft} ${jyear} ${cv}
    done
    exit
fi

if [ ${issh} -eq 1 ]; then
    CMD="python exec/ssh.py ${y1_clim} ${y2_clim}"
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
    export RUN=cp70
    export ORCA=ORCA1.L75
    export DIAG_D=/proj/bolinc/users/x_laubr/tmp/barakuda/ORCA1.L75_ece32b/ORCA1.L75-${RUN}
    export MM_FILE=/proj/bolinc/users/x_laubr/klaus/mesh_mask.nc
    export BM_FILE=/proj/bolinc/users/x_laubr/klaus/basin_mask.nc
    export NN_T="thetao"
    export NN_S="so"
    #
    cd ${DIAG_D}/
    python /home/x_laubr/DEV/barakuda/python/exec/plot_hovm_tz.py 1996 2000

    mv -f hov_*_ORCA1.L75-${RUN}*.png ${HERE}/
    #
fi








if [ ${iemp} -eq 1 ]; then

    export ORCA="ORCA1.L75"
    #export RUN="32bI"
    export RUN="cp00"
    export TSTAMP="1m"
    export DIAG_D="/proj/bolinc/users/x_laubr/tmp/barakuda/ORCA1.L75_ece32b/ORCA1.L75-${RUN}"
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





if [ ${icmip5} -eq 1 ]; then
    export RUN="SPIN"
    export CPREF="ORCA1-${RUN}_MM_"
    export ORCA=ORCA1.L42 ; # horizontal global configuration
    export NBL=42         ; # number of levels
    export STORE_DIR="/proj/bolinc/users/x_laubr"
    export TSTAMP="MM"
    export DIAG_D="."
    export MM_FILE="${STORE_DIR}/${ORCA}/${ORCA}-I/mesh_mask_ORCA1_ecearth2_42l_NoCaspian.nc4"
    export NN_SST="sosstsst"
    export NN_SSS="sosaline"
    export NN_SSH="sossheig"
    export NN_T="votemper"
    export NN_S="vosaline"

    export NN_TAUX="sozotaux"
    export NN_TAUY="sometauy"

    export FILE_DEF_BOXES="/home/x_laubr/DEV/barakuda/data/def_boxes_convection_ORCA1.txt"

    python /home/x_laubr/DEV/barakuda/python/exec/budget_rectangle_box.py 2250 100 uv

fi






if [ ${irnf} -eq 1 ]; then

    export ORCA="ORCA1.L75"
    export RUN="LB03"
    export DIAG_D="/proj/bolinc/users/x_laubr/tmp/barakuda/ORCA1.L75_ece32b/ORCA1.L75-${RUN}"
    export NN_RNF="runoffs"
    export MM_FILE="/proj/bolinc/users/x_laubr/tmp/barakuda/test/mesh_mask.nc"

    python exec/runoffs.py 1997 1999

fi




