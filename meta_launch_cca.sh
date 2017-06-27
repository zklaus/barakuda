#!/bin/bash -e

#---------------------------------------------------------------------------
# With PBSpro (qsub) we cannot pass argument to a script at the command line,
#  like this: 
#
#     qsub <qsub_options> script arg1 arg2
#
# Only this is posiible:
#
#     qsub <qsub_options> script
#
# So we rely on a template, which is parsed here to create a script that can
#  be submitted and has all the options we want.
#
#---------------------------------------------------------------------------
#    P. Le Sager, June/July 2017
#---------------------------------------------------------------------------

. src/bash/bash_functions.bash

usage() {
    barakuda_usage
    exit
}

# -------- Options
account=$ECE3_POSTPROC_ACCOUNT  # default account if exists
climato=0

while getopts C:R:i:e option ; do
    case $option in
        C) conf=${OPTARG} ;;
        R) exp=${OPTARG} ;;
        i) y1=${OPTARG} ;;
        e) y2=${OPTARG} ;;
        *) usage ;;
    esac
done
shift $((OPTIND-1))

[[ -n $y1 &&  -n $y2 ]] && cimato=1 # generate climato, and compare with it


# -------- Checks before submitting

[[ -z $conf ]] && usage
[[ -z $exp ]] && usage


# -- Scratch dir (location of submit script and its log)

OUT=$SCRATCH/tmp_barakuda
mkdir -p $OUT


# -------- Create diagnostics

cmd="./barakuda.sh -C ${conf} -R ${exp}"

tgt_script=$OUT/b_${exp}_diag.job

sed "s/<EXPID>/$exp/" < platform/cca.job.tmpl > $tgt_script

[[ -n $account ]] && \
    sed -i "s/<ACCOUNT>/$account/" $tgt_script || \
    sed -i "/<ACCOUNT>/ d" $tgt_script

sed -i "s|<DIAG>|diag|"  $tgt_script
sed -i "s|<OUT>|$OUT|"   $tgt_script
sed -i "s|<CMD>|${cmd}|" $tgt_script

diagid=$(qsub $tgt_script)
echo $diagid


# -------- Create climatologies

if (( $climato ))
then
    cmd="./build_clim.sh -C ${conf} -R ${exp} -i $y1 -e $y2"
    
    tgt_script=$OUT/b_${exp}_climato.job
    
    sed "s/<EXPID>/$exp/" < platform/cca.job.tmpl > $tgt_script
    
    [[ -n $account ]] && \
        sed -i "s/<ACCOUNT>/$account/" $tgt_script || \
        sed -i "/<ACCOUNT>/ d" $tgt_script
    
    sed -i "s|<DIAG>|clim|"  $tgt_script
    sed -i "s|<OUT>|$OUT|"   $tgt_script
    sed -i "s|<CMD>|${cmd}|" $tgt_script
    
    climid=$(qsub $tgt_script)
    echo $climid
fi


# -------- Create figures

if (( $climato ))
then
    cmd="./barakuda.sh -C ${conf} -R ${exp} -E"
else
    cmd="./barakuda.sh -C ${conf} -R ${exp} -e"
fi

tgt_script=$OUT/b_${exp}_fig.job

sed "s/<EXPID>/$exp/" < platform/cca.job.tmpl > $tgt_script

[[ -n $account ]] && \
    sed -i "s/<ACCOUNT>/$account/" $tgt_script || \
    sed -i "/<ACCOUNT>/ d" $tgt_script

sed -i "s|<DIAG>|fig|"  $tgt_script
sed -i "s|<OUT>|$OUT|"   $tgt_script
sed -i "s|<CMD>|${cmd}|" $tgt_script

if (( $climato ))
then
    qsub -W depend=afterok:$diagid:$climid  $tgt_script
else
    qsub -W depend=afterok:$diagid  $tgt_script
fi


# --------- info

qstat -wu $USER
