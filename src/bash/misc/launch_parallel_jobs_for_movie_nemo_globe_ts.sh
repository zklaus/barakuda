#!/bin/bash


DIR_NEMO_DATA=/scratch/Earth/lbrodeau/ORCA12-T1279/CHR0/extracted_nemo

#file=${DIR}/CHR0_1d_1990_${GRID}.nc4                                                                


#VAR=curl_ssu ; file=${DIR_NEMO_DATA}/CHR0_1d_1990_CURL.nc4
#VAR=sosstsst ; file=${DIR_NEMO_DATA}/${VAR}_CHR0_1990.nc
VAR=mod_ssu ; file=${DIR_NEMO_DATA}/CHR0_1d_1990_CURL.nc4


npj=20 ; istart=1 ; istop=99
#npj=20 ; istart=100 ; istop=199
#npj=20 ; istart=199 ; istop=300

#npj=8 ; istart=300 ; istop=364






icpt=1
jstrt=$((${istart}-1))

while [ $((${jstrt}+${npj})) -le $((${istop}+${npj})) ]; do
    jstop=$((${jstrt}+${npj}))


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
    jstrt=${jstop}
    icpt=$((${icpt}+1))
done
