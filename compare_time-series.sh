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

export script=compare_time-series
[ -z ${BARAKUDA_ROOT+x} ] && export BARAKUDA_ROOT=${PWD}

export FIG_FORMAT='png'
#export FIG_FORMAT='svg'


# Supported ORCA grids:
ORCA_LIST="ORCA025.L75 ORCA1.L75 ORCA1.L46 ORCA1.L42 ORCA2 ORCA2_L46"

# Display available configs:
list_conf="`\ls ${BARAKUDA_ROOT}/configs/config_*.sh | sed -e "s|${BARAKUDA_ROOT}/configs\/config_||g" -e s/'.sh'/''/g`"
# User configs, potentially in the directory from which barakuda.sh is called:
list_conf+=" `\ls ./config_*.sh 2>/dev/null | sed -e "s|.\/config_||g" -e s/'.sh'/''/g`"


# Important bash functions:
. ${BARAKUDA_ROOT}/src/bash/bash_functions.bash

usage()
{
    echo
    echo "USAGE: ${0} -C <config> -R <exp1,exp2,...,expN>  (options)"
    echo
    echo "     Available configs are:"
    for cc in ${list_conf}; do
        echo "         * ${cc}"
    done
    echo
    echo "   OPTIONS:"
    echo "      -y <YYYY> => force initial year to YYYY"
    echo
    echo "      -O <cnfg1,cnfg2,...,cnfgN> => list of ORCA configs corresponding to"
    echo "                                    <exp1,exp2,...,expN> in case they are"
    echo "                                    not all on the same ORCA grid"
    echo
#    echo "      -c <exp>  => 2D comparison diagnostics are performed against exp <exp>"
#    echo "                   instead of a climatology"
#    echo
    echo "      -f        => forces creation of diagnostic page eventhough treatment "
    echo "                    of output files is not finished"
    echo
    echo "      -e        => create the HTML diagnostics page on local or remote server"
    echo
    echo "      -h        => print this message"
    echo
    exit
}

CEXPS="" ; # list of experiments separated by ","
CORCS="" ; # corresponding list of ORCA configs in case CEXPS are not on the same ORCA grid!!!
YEAR0="" ; iforcey0=0

IPREPHTML=0
IFORCENEW=0

while getopts C:R:O:y:feh option ; do
    case $option in
        C) CONFIG=${OPTARG};;
        R) CEXPS=${OPTARG};;
        O) CORCS=${OPTARG};;
        y) YEAR0=${OPTARG} ; iforcey0=1 ;;
        f) IFORCENEW=1;;
        e) IPREPHTML=1;;
        h)  usage;;
        \?) usage ;;
    esac
done





if [ -z ${CONFIG} ] || [ -z ${CEXPS} ]; then usage ; exit ; fi



for og in ${ORCA_LIST}; do
    ca=""; ca=`echo ${CONFIG} | grep ${og}` ; if [ "${ca}" != "" ]; then ORCA=${og}; fi
done
if [ -z ${ORCA} ]; then echo "ORCA grid of config ${CONFIG} not supported yet"; exit; fi

echo

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
    echo "PROBLEM: cannot find file ${fconfig} !"; exit
fi
echo


LEXPS=`echo ${CEXPS} | sed -e s/'\,'/'\ '/g -e s/'\, '/'\ '/g`
echo; echo "Experiments to be treated: ${LEXPS}"
nbr=`echo ${LEXPS} | wc -w` ; echo "  => number of experiments to compare = ${nbr}"


if [ "${CORCS}" = "" ]; then
    if [ ! "${ORCA}" = "${CONF}" ]; then echo "ERROR: ORCA and CONF disagree! => ${ORCA} ${CONF}"; exit; fi
    export ORCA=${CONF}
else
    LORCS=`echo ${CORCS} | sed -e s/'\,'/'\ '/g -e s/'\, '/'\ '/g`
fi



# Should be set from bash_fucntions:
PYTH="${PYTHON_HOME}/bin/python -W ignore" ; # which Python installation to use
export PYTHONPATH=${PYTHON_HOME}/lib/python2.7/site-packages:${BARAKUDA_ROOT}/python/modules ; # PATH to python barakuda modules
PYBRKD_EXEC_PATH=${BARAKUDA_ROOT}/python/exec         ; # PATH to python barakuda executable
#-------------------


echo " NEMO grid = ${ORCA}";
echo " reading config into: ${fconfig}"
echo; echo

#NEXPS="${ORCA}-`echo ${LEXPS} | sed -e 's/\ /_/g'`"
NEXPS="`echo ${LEXPS} | sed -e 's/\ /_/g'`"
echo " Label to be used: ${NEXPS}" ; echo

BASE_NAME="comp_${NEXPS}"

DIAG_COMP_DIR=${DIAG_DIR}/comparisons/${BASE_NAME} ; rm -rf ${DIAG_COMP_DIR} ; mkdir -p ${DIAG_COMP_DIR}

YEAR_INI=4000
YEAR_END=0

# just that they become arrays...
VEXPS=( ${LEXPS} ) ;  VCONFEXPS=( ${LEXPS} ) ; VDIAGS=( ${LEXPS} ) ;

if [ "${CORCS}" = "" ]; then
    VORCS=( ${LEXPS} )
else
    VORCS=( ${LORCS} )
fi

jr=0

for exp in ${LEXPS}; do

    echo; echo " EXP ${exp} "
    EXP="${exp}"

    
    echo "   => ORCA = ${ORCA}"
    if [ "${CORCS}" = "" ]; then
        VORCS[${jr}]=${ORCA}
    fi


    CONFEXP=${VORCS[${jr}]}-${EXP}
    echo "   => CONFEXP = ${CONFEXP}"
    VCONFEXPS[${jr}]=${CONFEXP}


    DIAG_D="${DIAG_DIR}/${CONFEXP}"
    VDIAGS[${jr}]=${DIAG_D}
    echo "   => DIAG_D = ${DIAG_D} "; echo ; echo


    if [ ! -d ${DIAG_D} ]; then
        echo "PROBLEM: ${DIAG_D} does not exist!"
        echo "    =>  you must run barakuda for ${exp} prior to comparison!"
        exit
    fi


    # Guessing initial and last year:

    check_if_file ${DIAG_D}/first_year.info
    iy=`cat ${DIAG_D}/first_year.info`
    if [ ${iy} -lt ${YEAR_INI} ]; then export YEAR_INI=${iy}; fi

    check_if_file ${DIAG_D}/last_year_done.info
    iy=`cat ${DIAG_D}/last_year_done.info`
    if [ ${iy} -gt ${YEAR_END} ]; then export YEAR_END=${iy}; fi

    ((jr++))
done


echo
echo " Global YEAR_INI = ${YEAR_INI}"
echo " Global YEAR_END = ${YEAR_END}"
echo

#echo " VEXPS => ${VEXPS[*]} "
#echo " VORCS => ${VORCS[*]} "
#echo " VCONFEXPS => ${VCONFEXPS[*]} "
#echo " VDIAGS => ${VDIAGS[*]} "

export LIST_EXPS=${LEXPS}
export LIST_CONF=${VORCS[*]}

echo
echo "LIST_EXPS = ${LIST_EXPS}"
echo "LIST_CONF = ${LIST_CONF}"
echo



cd ${DIAG_COMP_DIR}/

${PYTH} ${PYBRKD_EXEC_PATH}/compare_time_series.py ${YEAR_INI} ${YEAR_END}


# Starting to configure HTML index file:
if [ "${EXTRA_CONF}" = "" ]; then echo "Problem, variable EXTRA_CONF is not set!" ; exit; fi
TITLE="Ocean diagnostics<br>Comparison of experiments: \"${VEXPS[*]}\"<br>Configuration: ${ORCA}_${EXTRA_CONF}"
if [ ${ece_exp} -gt 0 ]; then TITLE="${TITLE}<br>Atmospheric model: ${AGCM_INFO}"; fi

export CONFEXP="Comparison"

. ${BARAKUDA_ROOT}/src/bash/build_html.bash

parse_html ${BARAKUDA_ROOT}/src/html/conf_start.html index.html


list_figs=`\ls -v *.${FIG_FORMAT}`

for ff in ${list_figs}; do
    echo "<br><br><big> `echo ${ff} | sed -e s/.${FIG_FORMAT}//g -e s/_comparison//g` </big><br>" >> index.html
    echo "  <img style=\"border: 0px solid\" alt=\"\" src=\"${ff}\"> <br>" >> index.html

done

cat ${BARAKUDA_ROOT}/src/html/conf_end.html >> index.html ; # Closing HTML file...

cp ${BARAKUDA_ROOT}/src/html/conf_*.html  .
if [ ${ece_exp} -eq 0 ]; then
    cp ${BARAKUDA_ROOT}/src/html/logo.*g .
else
    cp ${BARAKUDA_ROOT}/src/html/logo_ece.svg ./logo.svg
    cp ${BARAKUDA_ROOT}/src/html/logo_ece.png ./logo.png
fi
echo; echo


if [ ${ihttp} -eq 0 ]; then

    echo "Diagnostic page installed in `pwd`"
    echo " => view this directory with a web browser (index.html)..."

else

    ssh ${RUSER}@${RHOST} "mkdir -p ${RWWWD}"

    if [ ${ihttp} -eq 1 ]; then

        echo "Preparing to export to remote host!"; echo

        cd ../

        tar cvf ${BASE_NAME}.tar ${BASE_NAME}
        scp ${BASE_NAME}.tar ${RUSER}@${RHOST}:${RWWWD}/
        ssh ${RUSER}@${RHOST} "cd ${RWWWD}/; rm -rf ${BASE_NAME}; tar xf ${BASE_NAME}.tar 2>/dev/null; rm -f ${BASE_NAME}.tar; chmod -R a+r ${BASE_NAME}"
        rm -f ${BASE_NAME}.tar

        echo; echo
        echo "Diagnostic page installed on remote host ${RHOST} in ${RWWWD}/${BASE_NAME}!"
        echo "( Also browsable on local host in `pwd` )"

    else

        echo "Error: \"ihttp\" is either 0 or 1 !"

    fi

fi

rm -rf ${COMP_DIR}

echo; echo
