#!/bin/bash

file=CHR0_1d_1990_CURL.nc4

#npj=20 ; istart=150 ; istop=365
npj=15 ; istart=350 ; istop=365

icpt=1
jstrt=${istart}

while [ $((${jstrt}+${npj})) -le $((${istop}+${npj})) ]; do
    jstop=$((${jstrt}+${npj}))


    echo
    CMD="~/DEV/barakuda/python/exec/movie_nemo_globe_ts.py ${file} socurl mesh_mask.nc ${jstrt} ${jstop}"

    echo ${CMD}
    
    csc=job_${icpt}.bash

    cat > ${csc} <<EOF
#!/bin/bash
#
#######
#SBATCH -w gustafson
#SBATCH -n 1
#SBATCH -J PROJPLOT
#SBATCH -t 5:50:00
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

