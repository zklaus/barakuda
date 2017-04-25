#!/usr/bin/env bash

function barakuda_first_last_years()
{
    # Autosub !!!
    
    export nbfpy=1
    export cmmdd1="0101"
    export cmmdd2="1231"
    export cyy2="<y>"
    if [ ! -z ${NCHNKS_Y} ]; then
        export nbfpy=${NCHNKS_Y}
        if [ ${NCHNKS_Y} -eq 2 ]; then
            echo "Chuncks of 6 months!"
            export cmmdd1="0701"
        elif [ ${NCHNKS_Y} -eq 4 ]; then
            echo "Chuncks of 3 months!"
            export cmmdd1="1101"
            export cmmdd2="0131"
            export cyy2="<y+1>"
        fi
    fi
    
    # Initial month (in tar file name)...
    if [ -z ${M_INI_EC} ]; then
        export M_INI_EC="01"
    fi

    cd ${NEMO_OUT_D}/
    
    YEAR_INI=${Y_INI_EC}
    
    if ${LFORCE_YINI}; then
        if [ ${YEAR0} -lt ${YEAR_INI_F} ]; then echo "ERROR: forced initial year is before first year!"; exit; fi
        export YEAR_INI=${YEAR0}
        echo " Initial year forced to ${YEAR_INI} !"
    fi
    export YEAR_INI_F=${YEAR_INI} ; # saving the year deduced from first file
    
    for ctest in fc0 fc00; do
        cc="MMO_${EXP}_${Y_INI_EC}${M_INI_EC}01_${ctest}_"
        ca=`\ls ${cc}* 2>/dev/null`
        if [ ! "${ca}" = "" ]; then cpr=${cc}; export cfc0="${ctest}"; fi
    done
    if [ "${cpr}" = "" ]; then
        echo "ERROR: can't find expected files in `pwd`!" ; echo "${cc}"; echo
        exit
    fi
    
    export YEAR_END=`\ls ${cpr}* | sed -e s/"${cpr}"/''/g | tail -1 | cut -c1-4`
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
    cpr="MMO_${EXP}_${Y_INI_EC}${M_INI_EC}01_${cfc0}_"
    
    cycm=`printf "%04d" $((${jyear}-1))`
    cycp=`printf "%04d" $((${jyear}+1))`

    echo " *** Going to extract year ${cyear} into:"
    echo "   ${TMP_DIR}"
    list=`\ls ${cpr}*${cyear}*.tar`
    nw=`echo ${list} | wc -w`
    #if [ ${nw} -gt ${nbfpy} ]; then echo "ERROR: more than ${nbfpy} tar file for year ${cyear}!"; exit; fi
    echo "  => ${list}"; echo
    rsync -avP ${list} ${TMP_DIR}/
    echo
    cd ${TMP_DIR}/
    for ff in ${list}; do tar -xvf ${ff}; done

    cprfx="${EXP}_${TSTAMP}_${cyear}"    
    if [ ${nbfpy} -gt 1 ] && [ ! "${M_INI_EC}" = "01" ]; then cprfx="${EXP}_${TSTAMP}"; fi

    rm -f ${list} ${EXP}_1d_* ${EXP}_*_diaptr.nc.gz ${EXP}_*_icemoa.nc.gz
    ls -l
    list=`\ls ${cprfx}*.nc.gz`
    for ff in ${list}; do
        echo "gunzip -f ${ff} &"
        gunzip -f ${ff} &
    done
    wait
    ls -l
    echo
    
    # Need to concatenate files if more than 1 NEMO file per year:
    if [ ${nbfpy} -gt 1 ]; then

        # If run not initialized the first of January, it's a fucking mess!
        if [ ! "${M_INI_EC}" = "01" ]; then
            nn=$((${M_INI_EC}-1)) ; nm=$((12/${nbfpy}))
            nr_next_y=$((${nm}-${nn})) ;    # in next year group file, read rec 1->${nr_next_y}
            nr_prev_y=$((${nr_next_y}+1)) ; # lolo:maybe wrong!!! in previous year group file, read rec nr_prev_y->12
            #
            for gt in ${NEMO_SAVED_FILES}; do
                echo "ncks -O -F -d time_counter,${nr_prev_y},${nm} ${cprfx}_${cycm}*_${cyear}*_${gt}.nc -o ${cprfx}_${cyear}0000_${gt}.nc"
                ncks -O -F -d time_counter,${nr_prev_y},${nm} ${cprfx}_${cycm}*_${cyear}*_${gt}.nc -o ${cprfx}_${cyear}0000_${gt}.nc
                rm -f ${cprfx}_${cycm}*_${cyear}*_${gt}.nc
                echo "ncks -O -F -d time_counter,1,${nr_next_y} ${cprfx}_${cyear}*_${cycp}*_${gt}.nc -o ${cprfx}_${cyear}9999_${gt}.nc"
                ncks -O -F -d time_counter,1,${nr_next_y} ${cprfx}_${cyear}*_${cycp}*_${gt}.nc -o ${cprfx}_${cyear}9999_${gt}.nc
                rm -f ${cprfx}_${cyear}*_${cycp}*_${gt}.nc
            done            
        fi
        
        echo; echo " *** Going to ncrcat !!!***"; echo
        for gt in ${NEMO_SAVED_FILES}; do
            \ls ${cprfx}*_${gt}.nc
            echo "ncrcat -O ${cprfx}*_${gt}.nc -o tmp_${gt}.nc &"
            ncrcat -O ${cprfx}*_${gt}.nc -o tmp_${gt}.nc &
            echo
        done
        wait
        for gt in ${NEMO_SAVED_FILES}; do
            rm -f ${cprfx}*_${gt}.nc
            mv -f tmp_${gt}.nc ${CRT1M}_${gt}.nc
        done
        rm -f ${cprfx}*.nc.gz
    fi
    
    # Testing if ALL required files are present now:
    for gt in ${NEMO_SAVED_FILES}; do
        ftt="./${CRT1M}_${gt}.nc" ;  check_if_file ${ftt}
    done
    echo; echo "All required files are in `pwd` for year ${cyear} !"; echo
}


function barakuda_check_year_is_complete()
{
    # Autosub !!!
    cycm=`printf "%04d" $((${jyear}-1))`
    cycp=`printf "%04d" $((${jyear}+1))`
    export i_get_file=1
    
    export TTAG=${cyear}0101_${cyear}1231 # calendar-related part of the file name    
    ttag_test=${cyear}${cmmdd1}-${cyy2}${cmmdd2}
    ttag_test=`echo ${ttag_test} | sed -e s/"<y>"/"${cyear}"/g | sed -e s/"<y+1>"/"${cycp}"/g` 

    # Testing if the current year-group has been done
    tdir=`echo ${NEMO_OUT_STRCT} | sed -e "s|<ORCA>|${ORCA}|g" -e "s|<EXP>|${EXP}|g"  -e "s|<Y_INI_EC>|${Y_INI_EC}|g"`
    cftar="${tdir}/MMO_${EXP}_${Y_INI_EC}${M_INI_EC}01_${cfc0}_${ttag_test}.tar"
    if ${lcontinue}; then
        if [ ! -f ${cftar} ]; then
            echo "Year ${cyear} is not completed yet:"; echo " => ${cftar} is missing"; echo
            export lcontinue=false
        fi
    fi
    if ${lcontinue}; then echo " *** Archive(s) for ${TTAG} is there!"; echo "  => (${cftar})"; fi
    echo
}
