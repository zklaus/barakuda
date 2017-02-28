#!/usr/bin/env bash

function build_index_html()
{
    echo; echo; echo
    echo "Creating HTML file!"

    ctl='<br><br><br><big><big>' ; ctr='</big></big><br><br>'
    csl='<br><br><big>' ; csr='</big><br>'
    spf='<br><br>'
    img_l='<img style="border: 0px solid" alt="" src="' ; img_r='"> <br><br>'

    cr="${CONFRUN}" ; ff="${FIG_FORM}"

    cd ${DIAG_D}/

    rm -f index.html
    
    if [ "${EXTRA_CONF}" = "" ]; then echo "Problem, variable EXTRA_CONF is not set!" ; exit; fi

    TITLE="Ocean diagnostics<br>Experiment: \"${RUN}\"<br>Configuration: ${ORCA}_${EXTRA_CONF}"
    if [ ${ece_run} -gt 0 ]; then TITLE="${TITLE}<br>Atmospheric model: ${AGCM_INFO}"; fi

    # Starting to configure HTML index file:
    sed -e "s|{TITLE}|${TITLE}|g" -e "s|{CONFRUN}|${CONFRUN}|g" -e "s|{DATE}|`date`|g" -e "s|{HOST}|${HOST}|g" \
        ${BARAKUDA_ROOT}/src/html/conf_start.html > index.html

    # Climato:
    if ${l_pclim}; then
        cat >> index.html <<EOF
    ${ctl} Diags from climatology (${CLIM_PER}) ${ctr}
    <big> <a href="./temp_sal/index.html"> Temperature and Salinity vs CLIM</a> </big>             ${spf}
        <big> <a href="./ssh/index.html">  Sea Surface Height </a> </big>                              ${spf}
        <big> <a href="./sea_ice/index.html">  Arctic and Antarctic sea-ice extent vs CLIM </a> </big> ${spf}
        <big> <a href="./mld/index.html">  Mixed Layer depth in relevent regions </a> </big>           ${spf}
        <big> <a href="./moc/index.html">  Meridional Overturning Circulation </a> </big>              ${spf}
EOF
        if [ ${i_do_sect} -eq 1 ]; then
            echo "<big> <a href='./temp_sal/index_sections.html'> Zonal/Meridional sections of T & S vs CLIM</a> </big>    ${spf}" >> index.html
        fi
        #
        if ${lcomp_to_run}; then
            cat >> index.html <<EOF
            ${ctl} Comparison with run ${RUNREF}, climatology (2004-2007) ${ctr}
            <big> <a href="./temp_sal/index_${RUNREF}.html"> Temperature and Salinity vs ${RUNREF}</a> </big>             ${spf}
            <!--        <big> <a href="./ssh/index_${RUNREF}.html">  Sea Surface Height </a> </big>                       ${spf}
            <big> <a href="./sea_ice/index_${RUNREF}.html">  Arctic and Antarctic sea-ice extent vs ${RUNREF} </a> </big> ${spf}
            <big> <a href="./mld/index_${RUNREF}.html">  Mixed Layer depth in relevent regions </a> </big>                ${spf}
            <big> <a href="./moc/index_${RUNREF}.html">  Meridional Overturning Circulation </a> </big>                   ${spf}
            -->
EOF
        fi
    fi

    if [ ${i_do_movi} -eq 1 ]; then
        cat >> index.html <<EOF
        ${ctl} Evolution of SST and SSS biases (w.r.t observations) ${ctr}
        ${img_l} dsst_${cr}.gif ${img_r}
        ${img_l} dsss_${cr}.gif ${img_r}
        ${img_l}  mld_${cr}.gif ${img_r}
EOF
    fi

    # AMOC page:
    cat >> index.html <<EOF
    ${ctl} Atlantic Meridional Overturning Circulation ${ctr}
    ${img_l} amoc_${cr}.${ff} ${img_r}
    ${img_l} amoc_${cr}_comp.${ff} ${img_r}
EOF

    # Temperature page
    cat >> index.html <<EOF
    ${ctl} Temperature time-series ${ctr}
    ${img_l} 3d_thetao_${cr}.${ff} ${img_r}
    ${img_l} mean_tos_${cr}.${ff} ${img_r}
    ${img_l} 3d_thetao_lev_${cr}.${ff} ${img_r}
    ${img_l} 3d_thetao_basins_${cr}.${ff} ${img_r}
    ${img_l} Nino34_${cr}.${ff} ${img_r}
EOF
    # AMO figure if here:
    famo="mean_SST_NAtl_${cr}.${ff}"
    if [ -f ${HTML_DIR}/${famo} ]; then echo "    ${img_l} ${famo} ${img_r}" >> index.html ; fi
    
    list_hov_figs=`\ls -v ${HTML_DIR}/hov_temperature_${cr}*.${ff}`
    if [ ! "${list_hov_figs}" = "" ]; then
        echo "    ${ctl} Time-depth evolution of temperature${ctr}" >> index.html
        for fhov in ${list_hov_figs}; do
            fgn=`basename ${fhov}`
            echo "    ${img_l} ${fgn} ${img_r}"  >> index.html
        done
    fi
    
    # Salinity page
    cat >> index.html <<EOF
    ${ctl} Salinity time-series ${ctr}
    ${img_l} 3d_so_${cr}.${ff} ${img_r}
    ${img_l} mean_sos_${cr}.${ff} ${img_r}
    ${img_l} 3d_so_lev_${cr}.${ff} ${img_r}
    ${img_l} 3d_so_basins_${cr}.${ff} ${img_r}
EOF
    list_hov_figs=`\ls -v ${HTML_DIR}/hov_salinity_${cr}*.${ff}`
    if [ ! "${list_hov_figs}" = "" ]; then
        echo "    ${ctl} Time-depth evolution of salinity${ctr}" >> index.html
        for fhov in ${list_hov_figs}; do
            fgn=`basename ${fhov}`
            echo "    ${img_l} ${fgn} ${img_r}"  >> index.html
        done
    fi
    

    # Surface heat flux diagnostics:
    if [ ! "${NN_QNET}" = "X" ]; then
        LIST_HF_FIG="htf_qnt htf_qsr \
        htf_qnt_NEMO_IFS htf_qnt_NEMO_IFS_annual \
        htf_qsr_NEMO_IFS htf_qsr_NEMO_IFS_annual"
        #
        cat >> index.html <<EOF
    ${ctl} Surface Heat flux time-series ${ctr}
EOF
        for fd in ${LIST_HF_FIG}; do
            fgn="mean_${fd}_${cr}.${ff}"; fgf="${HTML_DIR}/${fgn}"
            if [ -f ${fgf} ]; then
                echo "${img_l} ${fgn} ${img_r}" >> index.html
            fi
        done
    fi
    
    # Freshwater flux diagnostics:
    LIST_FW_FIG="zos fwf_fwf fwf_emp fwf_prc fwf_rnf fwf_clv fwf_rnf_clv \
        fwf_evp_NEMO_IFS fwf_evp_NEMO_IFS_annual \
        fwf_prc_NEMO_IFS fwf_prc_NEMO_IFS_annual \
        fwf_EmP_NEMO_IFS fwf_EmP_NEMO_IFS_annual \
        fwf_rnf_NEMO_IFS fwf_rnf_NEMO_IFS_annual \
        fwf_EmPmR_NEMO_IFS fwf_EmPmR_NEMO_IFS_annual \
        fwf_prc_IFS fwf_emp_ALL_IFS"
    #
    cat >> index.html <<EOF
    ${ctl} Surface Freshwater flux time-series ${ctr}
EOF
    for fd in ${LIST_FW_FIG}; do
        fgn="mean_${fd}_${cr}.${ff}"; fgf="${HTML_DIR}/${fgn}"
        if [ -f ${fgf} ]; then
            echo "${img_l} ${fgn} ${img_r}" >> index.html
        fi
    done



    # Sea-ice page
    if [ ${i_do_ice}  -gt 0 ]; then
        if [ ${i_do_movi} -eq 1 ]; then
            cat >> index.html <<EOF
            ${ctl} Evolution of Arctic/Antarctic concentration ${ctr}
            ${img_l} icen_${cr}.gif ${img_r}
            ${img_l} ices_${cr}.gif ${img_r}
EOF
        fi
        cat >> index.html <<EOF
        ${ctl} Arctic/Antarctic sea-ice time-series${ctr}
        ${img_l} seaice_extent_winter_${cr}.${ff} ${img_r}
        ${img_l} seaice_extent_summer_${cr}.${ff} ${img_r}
        ${img_l} seaice_volume_winter_${cr}.${ff} ${img_r}
        ${img_l} seaice_volume_summer_${cr}.${ff} ${img_r}
EOF
    fi

    if [ ${i_do_trsp} -gt 0 ]; then
        # Adding transport section part:
        echo "    ${ctl} Transport through sections ${ctr}" >> index.html
        list_section=`cat ${TRANSPORT_SECTION_FILE} | grep -v '^#' | grep '-'`
        for cs in ${list_section}; do
            echo ${cs}
            echo "    ${img_l} transport_vol_${cs}_${cr}.${ff} ${img_r}"  >> index.html
            echo "    ${img_l} transport_heat_${cs}_${cr}.${ff} ${img_r}" >> index.html
            echo "    <br>" >> index.html
        done
    fi

    # Checking if figures with time-series of MLD in specified boxes are here and adding them:
    if [ ${i_do_mean} -eq 1 ]; then
        list_mld_figs=`\ls -v ${HTML_DIR}/mean_mldr10_1_${cr}*.${ff}`
        if [ ! "${list_mld_figs}" = "" ]; then
            echo "    ${ctl} Horizontally-averaged Mixed-Layer Depth in different regions${ctr}" >> index.html
            for fmld in ${list_mld_figs}; do
                fgn=`basename ${fmld}`
                echo "    ${img_l} ${fgn} ${img_r}"  >> index.html
            done
        fi
    fi

    if [ ${i_do_sigt} -eq 1 ]; then
        # Adding transport by sigma class section part:
        echo "${ctl} Transport by sigma class at overflow sills${ctr}" >> index.html
        list_section=`cat ${DENSITY_SECTION_FILE} | grep -v '^#' | grep '_'`
        for cs in ${list_section}; do
            echo ${cs}
            echo "    ${img_l} transport_sigma_class_${cs}_${cr}.${ff} ${img_r}"  >> index.html
        done
        echo "    ${img_l} tr_sigma_gt278_${cr}.${ff} ${img_r}"  >> index.html
        echo "    ${spf}" >> index.html
    fi

    if [ ${i_do_mht} -eq 1 ]; then
        # Adding meridional heat transport:
        echo "${ctl} Meridional transports${ctr}"  >> index.html
        for coce in "GLO" "atl" "pac" "ind"; do
            echo "    ${img_l} MHT_${cr}_${coce}.${ff} ${img_r}"     >> index.html
            echo "    ${img_l} MST_${cr}_${coce}.${ff} ${img_r}" >> index.html
        done
        echo "    ${spf}" >> index.html
    fi

    cat ${BARAKUDA_ROOT}/src/html/conf_end.html >> index.html

    }




function build_sub_html()
{
    echo; echo; echo
    echo "Creating sub HTML files!"
        
    ctl='<br><br><br><big><big>' ; ctr='</big></big><br><br>'
    spf='<br><br>'
    img_l='<img style="border: 0px solid" alt="" src="' ; img_r='"> <br><br>'
    
    cr="${CONFRUN}" ; ff="${FIG_FORM}"

    cd ${DIAG_D}/
    
       # T, S, SSH and ice HTML page:
    for cdiag in ${DIRS_2_EXP}; do
        cat ${BARAKUDA_ROOT}/src/html/conf_start.html               > index.tmp
        cat ${BARAKUDA_ROOT}/src/html/${cdiag}/index_${cdiag}.html >> index.tmp
        cat ${BARAKUDA_ROOT}/src/html/conf_end.html                >> index.tmp
        sed -e "s|{TITLE}|${TITLE}|g" -e "s|{CONFRUN}|${cr}|g" \
            -e "s|{DATE}|`date`|g" -e "s|{HOST}|${HOST}|g" -e "s|{COMP2D}|CLIM|g" \
            index.tmp > ${cdiag}/index.html
        rm -f index.tmp
        cd ${cdiag}/ ; ln -sf ../logo.*g . ; cd ../
    done
    
    for var in "sst" "sss" "ts_100m" "ts_1000m" "ts_3000m"; do
        cat ${BARAKUDA_ROOT}/src/html/conf_start.html               > index.tmp
        cat ${BARAKUDA_ROOT}/src/html/temp_sal/${var}.html         >> index.tmp
        cat ${BARAKUDA_ROOT}/src/html/conf_end.html                >> index.tmp
        sed -e "s|{TITLE}|${TITLE}|g" -e "s|{CONFRUN}|${cr}|g" \
            -e "s|{DATE}|`date`|g" -e "s|{HOST}|${HOST}|g" -e "s|{COMP2D}|CLIM|g" \
            index.tmp > temp_sal/${var}_CLIM.html
    done
    
    # T&S sections:
    if [ ${i_do_sect} -eq 1 ]; then
        cat ${BARAKUDA_ROOT}/src/html/conf_start.html > index.tmp
        list_sec_figs=`\ls -v ${DIAG_D}/temp_sal/section_T_*_${cr}.${ff}`
        if [ ! "${list_sec_figs}" = "" ]; then
            echo "    ${ctl} Meridional and zonal cross-sections of T & S (specified in data/TS_sections.dat) ${ctr}" >> index.tmp
            for fsct in ${list_sec_figs}; do
                fgnt_n=`basename ${fsct}`                                      ; #figure for T NEMO
                fgnt_c=`echo ${fgnt_n} | sed -e s/"${cr}"/"CLIM"/g`            ; #figure for T CLIM
                fgns_n=`echo ${fgnt_n} | sed -e s/"section_T_"/"section_S_"/g` ; #figure for S NEMO
                fgns_c=`echo ${fgnt_c} | sed -e s/"section_T_"/"section_S_"/g` ; #figure for S CLIM
                cnames=`echo ${fgnt_n} | sed -e s/"section_T_"/""/g -e s/"_${cr}.${ff}"/""/g` ; # Name of section
                cnames=`echo ${cnames} | sed -e s/"_"/" "/g`
                echo "    ${csl} Cross-section '${cnames}': ${csr}"  >> index.tmp
                echo "    ${img_l} ${fgnt_n} ${img_r}"          >> index.tmp
                echo "    ${img_l} ${fgnt_c} ${img_r}"          >> index.tmp
                echo "    ${img_l} ${fgns_n} ${img_r}"          >> index.tmp
                echo "    ${img_l} ${fgns_c} ${img_r} <br><br><br>" >> index.tmp
            done
            cat ${BARAKUDA_ROOT}/src/html/conf_end.html          >> index.tmp
            sed -e "s|{TITLE}|${TITLE}|g" -e "s|{CONFRUN}|${cr}|g" \
                -e "s|{DATE}|`date`|g" -e "s|{HOST}|${HOST}|g" -e "s|{COMP2D}|CLIM|g" \
                index.tmp > temp_sal/index_sections.html
            rm -f index.tmp
        fi
    fi
    
    if ${lcomp_to_run}; then
        for cdiag in ${DIRS_2_EXP_RREF}; do
            cat ${BARAKUDA_ROOT}/src/html/conf_start.html               > index.tmp
            cat ${BARAKUDA_ROOT}/src/html/${cdiag}/index_${cdiag}.html >> index.tmp
            cat ${BARAKUDA_ROOT}/src/html/conf_end.html                >> index.tmp
            sed -e "s|{TITLE}|${TITLE}|g" -e "s|{CONFRUN}|${cr}|g" \
                -e "s|{DATE}|`date`|g" -e "s|{HOST}|${HOST}|g" -e "s|{COMP2D}|${RUNREF}|g" \
                index.tmp > ${cdiag}/index_${RUNREF}.html
            rm -f index.tmp
            cd ${cdiag}/ ; ln -sf ../logo.*g . ; cd ../
        done
        for var in "sst" "sss" "ts_100m" "ts_1000m" "ts_3000m"; do
            cat ${BARAKUDA_ROOT}/src/html/conf_start.html               > index.tmp
            cat ${BARAKUDA_ROOT}/src/html/temp_sal/${var}.html         >> index.tmp
            cat ${BARAKUDA_ROOT}/src/html/conf_end.html                >> index.tmp
            sed -e "s|{TITLE}|${TITLE}|g" -e "s|{CONFRUN}|${cr}|g" \
                -e "s|{DATE}|`date`|g" -e "s|{HOST}|${HOST}|g" -e "s|{COMP2D}|${RUNREF}|g" \
                index.tmp > temp_sal/${var}_${RUNREF}.html
        done
    fi
}



