#!/bin/bash -e

#---------------------------------------------------------------------------
# With PBSpro (qsub) you cannot pass arguments to a script at the command line,
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
# HOW-TO
#     you have to first "cd <location of barakuda package>/platform" and
#     then:
#
#       ./meta_launch_cca.sh -C ......
#
#---------------------------------------------------------------------------
# OPTIONAL SETTING
#
#   - you can overwrite the default hpc account for running barakuda with (in
#     your ~/.user_bashrc or ~/.bashrc):
#   
#       export ECE3_POSTPROC_ACCOUNT=<hpc account to use>
#
#     If not set, your default ECMWF account is used (i.e. 1st one in the list
#     you get with "account -l $USER" on ecgate). Note that you can also
#     specify an account at the command line when calling the script with the
#     -a option.
#
#---------------------------------------------------------------------------
#    P. Le Sager, June/July 2017
#---------------------------------------------------------------------------

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
    echo "   -a ACCOUNT    : specify a different special project for accounting (default: \$ECE3_POSTPROC_ACCOUNT)"
    echo "   -y YEAR       : year to start the processing (default: start year of the"
    echo "                        experiment, which is set in the config file!)"
    echo "   -i YEAR_START : first year to build climatolgy. Climatology is built if set."
    echo "   -e YEAR_END   :  last year to build climatolgy. Climatology is built if set."
    echo
    exit
}

# -------- Checks before submitting
[[ -z $ECE3_BARAKUDA_TOPDIR ]] && cd .. && ECE3_BARAKUDA_TOPDIR=$PWD
cd $ECE3_BARAKUDA_TOPDIR

# Available configs:
list_conf=$(\ls ${ECE3_BARAKUDA_TOPDIR}/configs/config_*.sh | sed -e "s|${ECE3_BARAKUDA_TOPDIR}/configs\/config_||g" -e s/'.sh'/''/g)

# User configs, potentially in the directory from which barakuda.sh is called:
list_conf+=" $(\ls ./config_*.sh 2>/dev/null | sed -e "s|.\/config_||g" -e s/'.sh'/''/g)"

# -------- Options
account=$ECE3_POSTPROC_ACCOUNT  # default account if exists
climato=0
opts=''

while getopts C:R:a:i:e:y: option ; do
    case $option in
        C) conf=${OPTARG} ;;
        R) exp=${OPTARG} ;;
        a) account=$OPTARG ;;
        i) y1=${OPTARG} ;;
        e) y2=${OPTARG} ;;
        y) opts="-y ${OPTARG}" ;;
        *) usage ;;
    esac
done
shift $((OPTIND-1))

[[ -n $y1 &&  -n $y2 ]] && climato=1 # generate climato, and compare with it

[[ -z $conf ]] && usage
[[ ! ${list_conf} =~ ${conf} ]] && { echo ; echo "!!!! UNKNOWN CONF: $conf !!!!!"; usage;}
[[ -z $exp ]] && usage


# -- Scratch dir (location of submit script and its log)

OUT=$SCRATCH/tmp_barakuda
mkdir -p $OUT


# -------- Create diagnostics

cmd="./barakuda.sh ${opts} -C ${conf} -R ${exp}"

tgt_script=$OUT/b_${exp}_diag.job

sed "s/<EXPID>/$exp/" < platform/cca.job.tmpl > $tgt_script

[[ -n $account ]] && \
    sed -i "s/<ACCOUNT>/$account/" $tgt_script || \
    sed -i "/<ACCOUNT>/ d" $tgt_script

sed -i "s|<BARAKUDA_TOPDIR>|${ECE3_BARAKUDA_TOPDIR}|"  $tgt_script
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
    
    sed -i "s|<BARAKUDA_TOPDIR>|${ECE3_BARAKUDA_TOPDIR}|"  $tgt_script
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

sed -i "s|<BARAKUDA_TOPDIR>|${ECE3_BARAKUDA_TOPDIR}|"  $tgt_script
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
