#!/usr/bin/env bash

function barakuda_import_files()
{

    echo " SPECIFIC 1y 3D !!! into `pwd` !"

    CROUT=${CPREF}${TTAG}
    echo " CROUT => ${CROUT}"

    CRTY=`echo ${CROUT} | sed -e s/"_1m_"/"_${ANNUAL_3D}_"/g`

    echo " CRTY => ${CRTY}"

    # Import command:
    CIMP="rsync -L"
    if [ "${CONF}" = "ORCA025.L75" ]; then CIMP="ln -sf" ; fi

    # On what file type to test file presence:
    cgrid_test_1m=`echo ${NEMO_SAVED_FILES} | cut -d ' ' -f2`
    cgrid_test_1y=`echo ${NEMO_SAVED_FILES_3D} | cut -d ' ' -f2`

    echo " *** testing on files \"${cgrid_test_1m}\" and \"${cgrid_test_1y}\" !"; echo

#AB-----------------
    MYNEMO_OUT_D=`echo ${NEMO_OUT_STRCT} | sed -e "s|<EXP>|${EXP}|g" -e "s|????|${cyear}|g"`
#AB-----------------

    l_happy=false
    while ! ${l_happy} ; do
        if [ ${IFREQ_SAV_YEARS} -eq 1 ]; then l_happy=true; fi
        rm -f *.tmp
        if [ ${i_get_file} -eq 1 ]; then
            echo " => gonna get ${CROUT}_* files..."


            # Importing required "1m" files to tmp dir and unzipping:
            for gt in ${NEMO_SAVED_FILES}; do
                f2i=${CROUT}_${gt}.nc ;   sgz=""
                for ca in ".gz" "4"; do
#AB                    if [ -f ${NEMO_OUT_D}/${cpf}${f2i}${ca} ]; then sgz="${ca}"; fi
                    if [ -f ${MYNEMO_OUT_D}/${cpf}${f2i}${ca} ]; then sgz="${ca}"; fi
                done
#AB                check_if_file ${NEMO_OUT_D}/${cpf}${f2i}${sgz}
                check_if_file ${MYNEMO_OUT_D}/${cpf}${f2i}${sgz}
                if [ ! -f ./${f2i} ]; then
                    echo "Importing ${f2i}${sgz} ..."
#AB                    echo "${CIMP} ${NEMO_OUT_D}/${cpf}${f2i}${sgz} `pwd`/"
#AB                    ${CIMP} ${NEMO_OUT_D}/${cpf}${f2i}${sgz} ./
                    echo "${CIMP} ${MYNEMO_OUT_D}/${cpf}${f2i}${sgz} `pwd`/"
                    ${CIMP} ${MYNEMO_OUT_D}/${cpf}${f2i}${sgz} ./
                    if [ "${sgz}" = ".gz" ]; then gunzip -f ./${f2i}.gz ; fi
                    if [ "${sgz}" = "4"   ]; then
                        echo "mv -f ./${f2i}4 ./${f2i}"
                        mv -f ./${f2i}4 ./${f2i}
                    fi
                    check_if_file ${f2i}
                    echo " ... done!"; echo
                else
                    echo " ${f2i}${sgz} was already in `pwd`"
                fi
            done


            # Importing required "1y" files to tmp dir and unzipping:
            for gt in ${NEMO_SAVED_FILES_3D}; do
                f2i=${CRTY}_${gt}.nc ;   sgz=""
                for ca in ".gz" "4"; do
#AB                    if [ -f ${NEMO_OUT_D}/${cpf}${f2i}${ca} ]; then sgz="${ca}"; fi
                    if [ -f ${MYNEMO_OUT_D}/${cpf}${f2i}${ca} ]; then sgz="${ca}"; fi
                done
#AB                check_if_file ${NEMO_OUT_D}/${cpf}${f2i}${sgz}
                check_if_file ${MYNEMO_OUT_D}/${cpf}${f2i}${sgz}
                if [ ! -f ./${f2i} ]; then
                    echo "Importing ${f2i}${sgz} ..."
#AB                    echo "${CIMP} ${NEMO_OUT_D}/${cpf}${f2i}${sgz} `pwd`/"
#AB                    ${CIMP} ${NEMO_OUT_D}/${cpf}${f2i}${sgz} ./
                    echo "${CIMP} ${MYNEMO_OUT_D}/${cpf}${f2i}${sgz} `pwd`/"
                    ${CIMP} ${MYNEMO_OUT_D}/${cpf}${f2i}${sgz} ./
                    if [ "${sgz}" = ".gz" ]; then gunzip -f ./${f2i}.gz ; fi
                    if [ "${sgz}" = "4"   ]; then
                        echo "mv -f ./${f2i}4 ./${f2i}"
                        mv -f ./${f2i}4 ./${f2i}
                    fi
                    check_if_file ${f2i}
                    echo " ... done!"; echo
                else
                    echo " ${f2i}${sgz} was already in `pwd`"
                fi
            done

        fi

    done
    echo

    # Testing if ALL required files are present now:
    for gt in ${NEMO_SAVED_FILES}; do
        ftt="./${CROUT}_${gt}.nc" ;  check_if_file ${ftt}
    done
    for gt in ${NEMO_SAVED_FILES_3D}; do
        ftt="./${CRTY}_${gt}.nc" ;  check_if_file ${ftt}
    done
    echo; echo "All required files are in `pwd` for year ${cyear} !"; echo

}
