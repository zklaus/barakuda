#!/usr/bin/env bash

function barakuda_first_last_years()
{
    # Autosub !!!
    cd ${NEMO_OUT_D}/
    cpr="MMO_${EXP}_${Y_INI_EC}0101_fc00_"
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
    export NEMO_OUT_D=`echo ${NEMO_OUT_STRCT} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<EXP>|${EXP}|g" -e "s|<Y_INI_EC>|${Y_INI_EC}|"g`
    cd ${NEMO_OUT_D}/
    cpr="MMO_${EXP}_${Y_INI_EC}0101_fc00_"
    #
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
    rm -f ${list} ${EXP}_1d_* ${EXP}_*_diaptr.nc.gz
    ls -l
    gunzip -f ${EXP}_${TSTAMP}_${cyear}*.nc.gz
    ls -l
    echo
    #
    # Testing if ALL required files are present now:
    for gt in ${NEMO_SAVED_FILES}; do
        ftt="./${CRT1M}_${gt}.nc" ;  check_if_file ${ftt}
    done
    echo; echo "All required files are in `pwd` for year ${cyear} !"; echo
    #
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
    tdir=`echo ${NEMO_OUT_STRCT} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<EXP>|${EXP}|g"  -e "s|<Y_INI_EC>|${Y_INI_EC}|g"`
    cftar="${tdir}/MMO_${EXP}_${Y_INI_EC}0101_fc00_${cy1}0101-${cy1}1231.tar"
    if ${lcontinue}; then
        if [ ! -f ${cftar} ]; then
            echo "Year ${cy1} is not completed yet:"; echo " => ${cftar} is missing"; echo
            export lcontinue=false
        fi
    fi
    if ${lcontinue}; then echo " *** Archive for ${TTAG} is there!"; echo "  => (${cftar})"; fi
    echo
}
