#!/bin/bash


VAR=socurl ; GT=CURL; DIR_NEMO_DATA="."

FRQ="1d"




#npj=20 ; istart=150 ; istop=365
#npj=73 ; istart=0 ; istop=364
#npj=2 ; istart=69 ; istop=73
#npj=2 ; istart=137 ; istop=146
#npj=2 ; istart=206 ; istop=219

npj=5 ; istart=290 ; istop=365



file=${DIR_NEMO_DATA}/CHR0_${FRQ}_1990_${GT}.nc4


icpt=1
jstrt=${istart}

while [ $((${jstrt}+${npj})) -le $((${istop}+${npj})) ]; do
    jstop=$((${jstrt}+${npj}-1))


    echo
    CMD="~/DEV/barakuda/python/exec/movie_nemo_globe_ts.py ${file} ${VAR} mesh_mask.nc ${jstrt} ${jstop}"

    echo ${CMD}
    
    csc=job_${icpt}.sub

    cat > ${csc} <<EOF
#!/bin/bash
#
#######
#SBATCH -w gustafson
#SBATCH -n 2
#SBATCH -J PROJPLOT
#SBATCH -t 23:50:00
#SBATCH -o out_plot_globe_%J.out
#SBATCH -e err_plot_globe%J.err
########
export DIR_NCVIEW_CMAP=/home/Earth/lbrodeau/DEV/barakuda/src/ncview_colormaps
export PYTHONWARNINGS="ignore::DeprecationWarning"
${CMD}
EOF
    
    echo "sbatch ./${csc}"
    sbatch ./${csc}
    sleep 1    
    jstrt=$((${jstop}+1))
    icpt=$((${icpt}+1))
done
