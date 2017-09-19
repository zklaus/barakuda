#!/bin/bash
#
# based on NSC script of K. Wyser, July 2016
#

set -ex

usage() {
    echo
    echo "USAGE: ${0} -C <config> -R <experiment>  (options)"
    echo
    echo "     Available configs are:"
    for cc in ${list_conf}; do
        echo "         * ${cc}"
    done
    echo
    echo "Options are:"
    echo
    echo "   -i YEAR_START : first year to build climatolgy. Climatology is built if set."
    echo "   -e YEAR_END   :  last year to build climatolgy. Climatology is built if set."
    echo
    exit
}

climato=0
while getopts C:R:i:e: option ; do
    case $option in
        C) conf=${OPTARG} ;;
        R) expname=${OPTARG} ;;
        i) y1=${OPTARG} ;;
        e) y2=${OPTARG} ;;
        *) usage ;;
    esac
done
shift $((OPTIND-1))

[[ -n $y1 &&  -n $y2 ]] && climato=1 # generate climato, and compare with it

# -------- Checks before submitting
cd ..

# Available configs
list_conf=$(\ls ./configs/config_*.sh | sed -e "s|./configs\/config_||g" -e s/'.sh'/''/g)

# User configs, potentially in the directory from which barakuda.sh is called
list_conf+=" $(\ls ./config_*.sh 2>/dev/null | sed -e "s|.\/config_||g" -e s/'.sh'/''/g)"

[[ -z $conf ]] && usage
[[ -z $expname ]] && usage


job_name=b-$expname
launch_cmd='sbatch -n 1 -c 12'

# -------- Create diagnostics

ll=$( $launch_cmd -J $job_name-1 barakuda.sh -C $conf -R $expname )
echo $ll
job1_id=$(echo $ll | awk '{print $4}')

# -------- Create climatologies

if (( $climato ))
then
    ll=$( $launch_cmd -J $job_name-2 build_clim.sh -C $conf -R $expname -i $y1 -e $y2)
    echo $ll
    job2_id=$(echo $ll | awk '{print $4}')
fi

# -------- Create figures

if (( $climato ))
then
    $launch_cmd -d afterok:$job1_id:$job2_id -J $job_name-3 barakuda.sh -C $conf -R $expname -E
else
    $launch_cmd -d afterok:$job1_id -J $job_name-3 barakuda.sh -C $conf -R $expname -e
fi

