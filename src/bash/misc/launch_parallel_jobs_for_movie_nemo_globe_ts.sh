#!/bin/bash

file=CHR0_1d_1990_CURL.nc4 ; var=socurl

#npj=20 ; istart=150 ; istop=365
#npj=73 ; istart=0 ; istop=364
#npj=2 ; istart=69 ; istop=73
#npj=2 ; istart=137 ; istop=146
npj=2 ; istart=206 ; istop=219

icpt=1
jstrt=${istart}

while [ $((${jstrt}+${npj})) -le $((${istop}+${npj})) ]; do
    jstop=$((${jstrt}+${npj}))


    echo
    CMD="~/DEV/barakuda/python/exec/movie_nemo_globe_ts.py ${file} ${var} mesh_mask.nc ${jstrt} ${jstop}"

    echo ${CMD}
    
    csc=job_${icpt}.sub

    cat > ${csc} <<EOF
#!/bin/bash
#
#######
#SBATCH -w gustafson
#SBATCH -n 2
#SBATCH -J PROJPLOT
#SBATCH -t 11:50:00
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
    jstrt=$((${jstop}))
    icpt=$((${icpt}+1))
done
