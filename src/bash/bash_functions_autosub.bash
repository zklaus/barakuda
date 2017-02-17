#!/usr/bin/env bash


function barakuda_setup()
{
    script=`basename $0 | sed -e s/'.sh'/''/g`
    echo
    if [ ! "${ORCA}" = "${CONF}" ]; then echo "ERROR: ORCA and CONF disagree! => ${ORCA} ${CONF}"; exit; fi
    export ORCA=${CONF}
    echo

    if [ -z ${PYTHON_HOME} ]; then echo "ERROR: PYTHON_HOME is not set! => add it to config file"; exit; fi
    export PYTH="${PYTHON_HOME}/bin/python -W ignore" ; # which Python installation to use
    export PYTHONPATH=${PYTHON_HOME}/lib/python2.7/site-packages:${BARAKUDA_ROOT}/python/modules ; # PATH to python barakuda modules
    export PYBRKD_EXEC_PATH=${BARAKUDA_ROOT}/python/exec         ; # PATH to python barakuda executable

    echo " PYTHON_HOME => "${PYTHON_HOME} ; echo
    echo "  TRANSPORT_SECTION_FILE => ${TRANSPORT_SECTION_FILE} !" ; echo

    if ${l_clim_diag} ; then
        echo
        echo " Files containing climatologies to be used:"
        echo " T 3D => ${F_T_CLIM_3D_12} => ${NN_T_CLIM}"
        echo " S 3D => ${F_S_CLIM_3D_12} => ${NN_S_CLIM}"
        echo " SST  => ${SST_CLIM_12} => ${NN_SST_CLIM}"
        echo
        for ff in ${F_T_CLIM_3D_12} ${F_S_CLIM_3D_12} ${SST_CLIM_12}; do
            if [ ! -f ${F_T_CLIM_3D_12} ]; then echo "ERROR: ${ff} is missing!"; exit; fi
        done
    fi

    # Names for temperature, salinity, u- and v-current...
    if [ -z ${NN_T} ] || [ -z ${NN_S} ] || [ -z ${NN_U} ] || [ -z ${NN_V} ]; then
        echo "NN_T, NN_S, NN_U and NN_V are NOT given a value into"
        echo " in ${fconfig} "
        echo "  => using default names: thetao, so, uo, vo" ; echo
        NN_T="thetao"; NN_S="so"; NN_U="uo"; NN_V="vo"
    fi

    echo ; echo " *** NN_T=${NN_T}, NN_S=${NN_S}, NN_U=${NN_U} and NN_V=${NN_V} "; echo

    # Checking what files we have / plan to use:
    if [ -z "${NEMO_SAVED_FILES}" ]; then
        echo "Please specify which NEMO files are saved (file suffixes, grid_T, ..., icemod) ?"
        echo " => set the variable NEMO_SAVED_FILES in your config_${CONFIG}.sh file!"; exit
    fi
    echo; echo "File types to import (NEMO_SAVED_FILES) : ${NEMO_SAVED_FILES}"; echo; echo

    # Need to be consistent with the netcdf installation upon which cdftools_light was compiled:
    ff="cdftools_light/make.macro"
    if [ ! -f ${ff} ]; then echo "PROBLEM: cannot find ${ff} (needed to get NCDF_DIR)!"; exit; fi
    export NCDF_DIR=`cat ${ff} | grep ^NCDF_DIR | cut -d = -f2 | sed -e s/' '//g`
    echo ; echo "NCDF_DIR = ${NCDF_DIR}"; echo
    export LD_LIBRARY_PATH=${NCDF_DIR}/lib:${LD_LIBRARY_PATH}

    # Exporting some variables needed by the python scripts:
    export RUN=${RUN}
    export CONFRUN=${ORCA}-${RUN}
    export DIAG_D=${DIAG_DIR}/${CONFRUN}
    export CLIM_DIR=${DIAG_D}/clim

    # We need a scratch/temporary directory to copy these files to and gunzip them:
    if [ ! -d "${SCRATCH}" ]; then
        echo " ERROR: scratch directory ${SCRATCH} not there!!!"
        exit
    fi
    export TMP_DIR=${SCRATCH}/barakuda_${RUN}_tmp
    echo " IMPORTANT the TMP_DIR work directory is set to:" ; echo " ${TMP_DIR}"; echo ; sleep 2
    
    rm -rf ${TMP_DIR}
    mkdir -p ${DIAG_D} ${TMP_DIR}
    
    export NEMO_OUT_D=`echo ${NEMO_OUT_STRCT} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<RUN>|${RUN}|g"`
    if [ ! -d ${NEMO_OUT_D} ]; then echo "Unfortunately we could not find ${NEMO_OUT_D}"; exit; fi

    echo; echo " * Config to be used: ${CONFIG} => ORCA grid is ${ORCA}"
    echo " * Run is ${RUN} "; echo " * Files are stored into ${NEMO_OUT_D}"; echo; sleep 2

    export CPREF=`echo ${NEMO_FILE_PREFIX} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<RUN>|${RUN}|g" -e "s|<TSTAMP>|${TSTAMP}|g"`
    echo " NEMO files prefix = ${CPREF} "

    # only neede for barakuda.sh :
    if [ "${script}" = "barakuda" ]; then
        # Testing if NCO is installed:
        which ncks 1>out 2>/dev/null; ca=`cat out`; rm -f out
        if [ "${ca}" = "" ]; then echo "Install NCO!!!"; echo; exit; fi

        # Not fully supported yet:
        ca=" => diagnostic totally beta and not fully supported yet!"
        if [ ${i_do_amo}  -gt 0 ]; then echo " *** i_do_amo  ${ca}"; exit; fi
        if [ ${i_do_icet} -gt 0 ]; then echo " *** i_do_icet ${ca}"; exit; fi

        if [ ${ISTAGE} -eq 1 ]; then
            for ex in ${L_EXEC}; do check_if_file cdftools_light/bin/${ex} "Compile CDFTOOLS executables!"; done
        fi
    fi

    if [ -z ${NCDF_DIR} ]; then
        if [ ! -z ${NETCDF_DIR} ]; then
            export NCDF_DIR=${NETCDF_DIR}
        elif [ ! -z ${NETCDF_HOME} ]; then
            export NCDF_DIR=${NETCDF_HOME}
        else
            echo "ERROR: NCDF_DIR could not be determined, please specify it in your config file!"
            exit
        fi
    fi
    echo
}


function barakuda_first_last_years()
{
    # Autosub !!!
    cd ${NEMO_OUT_D}/
    cpr="MMO_${RUN}_${Y_INI_EC}0101_fc00_"
    YEAR_INI=${Y_INI_EC}
    #
    if ${LFORCE_YINI}; then
        if [ ${YEAR0} -lt ${YEAR_INI_F} ]; then echo "ERROR: forced initial year is before first year!"; exit; fi
        export YEAR_INI=${YEAR0}
        echo " Initial year forced to ${YEAR_INI} !"
    fi
    export YEAR_INI_F=${YEAR_INI} ; # saving the year deduced from first file
    #
    export YEAR_END=`\ls ${cpr}* | sed -e s/"${cpr}"/''/g | tail -1 | cut -c1-4`
    echo "LOLO YEAR_END => ${YEAR_END}"
    echo ${YEAR_END} |  grep "[^0-9]" >/dev/null; # Checking if it's an integer
    if [ ! "$?" -eq 1 ]; then
        echo "ERROR: it was imposible to guess the year coresponding to the last saved year!"
        echo "       => check your NEMO output directory and file naming..."; exit
    fi
    echo
    echo " *** Initial year set to ${YEAR_INI}"
    echo " ***   Last  year set to ${YEAR_END}"
    echo
}

function barakuda_import_files()
{
    # Autosub !!!
    echo
    echo "Inside barakuda_extract_autosub => year = ${cyear}"
    export NEMO_OUT_D=`echo ${NEMO_OUT_STRCT} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<RUN>|${RUN}|g"`
    cd ${NEMO_OUT_D}/
    cpr="MMO_${RUN}_${Y_INI_EC}0101_fc00_"
    #export TMP_DIR=${TMP_DIR}/NEMO ; mkdir -p ${TMP_DIR}
    echo " *** Going to extract year ${cyear} into:"
    echo "   ${TMP_DIR}"
    list=`\ls ${cpr}${cyear}*.tar`
    nw=`echo ${list} | wc -w`
    if [ ${nw} -gt 1 ]; then echo "ERROR: more than 1 tar file for year ${cyear}!"; exit; fi
    echo "  => ${list}"; echo
    rsync -avP ${list} ${TMP_DIR}/
    echo
    cd ${TMP_DIR}/
    tar -xvf ${list}
    rm -f ${list} ${RUN}_1d_* ${RUN}_*_diaptr.nc.gz
    ls -l
    gunzip -f ${RUN}_${TSTAMP}_${cyear}*.nc.gz
    ls -l
    echo
    #export NEMO_OUT_D=${TMP_NEMO_DIR}
}


function barakuda_check_year_is_complete()
{
    # Autosub !!!
    jy1=${jyear} ; jy2=$((${jyear}+${IFREQ_SAV_YEARS}-1))
    cy1=`printf "%04d" ${jy1}` ; cy2=`printf "%04d" ${jy2}`
    cy1m=`printf "%04d" $((${jy1}-${IFREQ_SAV_YEARS}))` ; cy2m=`printf "%04d" $((${jy2}-${IFREQ_SAV_YEARS}))`
    export i_get_file=1
    echo " *** (${cyear} => from ${cy1} to ${cy2})"
    export TTAG=${cy1}0101_${cy2}1231 # calendar-related part of the file name
    #
    # Testing if the current year-group has been done
    tdir=`echo ${NEMO_OUT_STRCT} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<RUN>|${RUN}|g"`
    cftar="${tdir}/MMO_${RUN}_${Y_INI_EC}0101_fc00_${cy1}0101-${cy1}1231.tar"
    if ${lcontinue}; then
        if [ ! -f ${cftar} ]; then
            echo "Year ${cy1} is not completed yet:"; echo " => ${cftar} is missing"; echo
            export lcontinue=false
        fi
    fi
    if ${lcontinue}; then echo " *** Archive for ${TTAG} is there!"; echo "  => (${cftar})"; fi
    echo
}
