#!/bin/bash

#==========================================================
#
#         Configuration file for
#
# OCEAN MONITORING for NEMO v2.? and EC-Earth 2.3 (CMIP5)
#
#            Machine: triolith
#
#        L. Brodeau, November 2016
#
#===========================================================

export CONF=ORCA1.L42 ; # horizontal global configuration
export NBL=42         ; # number of levels

export HOST=TRIOLITH ; # this has no importance at all, it will just become an "info" on the web-page!
export JTITLE="LIM2, NEMO 2.X (EC-Earth 2_CMIP5)" ;   #  // same here ...

# File system / path on which most netcdf data will be read:
export STORE_DIR="/data1/laurent"

# Path to directory containing some 2D and 3D climatologies on the relevant ORCA grid:
export CONF_INI_DIR="${STORE_DIR}/${CONF}/${CONF}-I"

# In what directory of the local machine to save the diagnostics:
export DIAG_DIR="/home/laurent/tmp/barakuda"

# --- only for problematic hosts ----
#module add NCO/4.2.3    
#module add PYTHON/2.7.3
# -----------------------------------

export PYTHON_HOME="/opt/Canopy_64bit/User" ; # HOME to python distribution with matplotlib and basemap

# Is it an ec-earth run?
export ece_run=0 ; # 0 => not an EC-Earth run, it's a "pure" ocean-only NEMO run done from traditional NEMO setup
#                  # 1 => it's an OCEAN-ONLY EC-Earth run done from a EC-Earth setup
#                  # 2 => it's a  COUPLED  EC-Earth run
#                  #      Both 1 and 2 imply that NEMO files are stored in something like
#                  #       ${STORE_DIR}/<RUN>/output/nemo/<YYY>
#                  #       where YYY starts from '001' to
#                  #   If you select '2', make sure 'cdo' is available and working!!!
#
export Y_INI_EC=1990 ;    # initial year if ec-earth run...

# List of suffixed of files that have been saved by NEMO and that are needed for the diags:
export NEMO_SAVED_FILES="grid_T grid_U grid_V icemod"

# Directory structure in which to find NEMO output file (use <ORCA> and <RUN>):
export NEMO_OUT_STRCT="${STORE_DIR}/CMIP5/<RUN>-SAVED/NEMO"

export TSTAMP="MM"   ; # output time-frequency stamp as in NEMO output files...

# How does the nemo files prefix looks like
# Everything before "<year_related_info>_grid_<X>" or "<year_related_info>_icemod"
# use <ORCA>, <RUN> and <TSTAMP>=>  Ex: export NEMO_FILE_PREFIX="<ORCA>-<RUN>_<TSTAMP>_"
export NEMO_FILE_PREFIX="ORCA1-<RUN>_<TSTAMP>_"
# => should get rid of TSTAMP actually...

# Temporary file system (scratch) on which to perform the job you can use <JOB_ID> if scracth depends on JOB ID:
export SCRATCH="/scratch/local/<JOB_ID>"


####### NEMO => what fields in what files ??? ############
#       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   => depends on the XIOS *.xml setup you used...
#   => always specify a string for the NN_* variables
#      even when missing from your files (ex: NN_MLD="xx")
#
# State variables and others in grid_T files:
export NN_SST="sosstsst"
export NN_SSS="sosaline"
export NN_SSH="sossheig"
export NN_T="votemper"
export NN_S="vosaline"
export NN_MLD="somxl010"
#
# State variables and others in grid_U files:
export NN_U="vozocrtx"
export NN_TAUX="sozotaux"
export NN_U_EIV="0" ; # 0 => ignore
# State variables and others in grid_V files:
export NN_V="vomecrty"
export NN_TAUY="sometauy"
export NN_V_EIV="0" ; # 0 => ignore
#
# Sea-ice fields:
export FILE_ICE_SUFFIX="icemod" ; # in what file type extension to find ice fields
export NN_ICEF="ileadfra" ; # name of ice fraction in "FILE_ICE_SUFFIX" file...
export NN_ICET="iicethic" ; # ice thickness but 'sit' is only in icemod file !!!
#
# Surface fluxes:
export FILE_FLX_SUFFIX="grid_T" ; # in what file type extension to find surface fluxes
export NN_FWF="sowaflup"        ; # name of net freshwater flux (E-P-R) in "FILE_FLX_SUFFIX" file...
export NN_EMP="XXX"             ; # name of E-P in "FILE_FLX_SUFFIX" file...
export NN_P="XXX"               ; # name of total precipitation (solid+liquid) in "FILE_FLX_SUFFIX" file...
export NN_RNF="sorunoff"        ; # name of continental runoffs in "FILE_FLX_SUFFIX" file...
export NN_CLV="XXX"             ; # calving from icebergs in "FILE_FLX_SUFFIX" file...
export NN_E="XXX"               ; # name of total evaporation in "FILE_FLX_SUFFIX" file...
export NN_QNET="sohefldo"       ; # name of total net surface heat flux in "FILE_FLX_SUFFIX" file...
export NN_QSOL="soshfldo"       ; # name of net surface solar flux in "FILE_FLX_SUFFIX" file...
#
################################################################################################

# Land-sea mask and basins files:
export MM_FILE="${CONF_INI_DIR}/mesh_mask_ORCA1_ecearth2_42l_NoCaspian.nc4"
export BM_FILE="${CONF_INI_DIR}/basin_mask_ORCA1_ecearth2_42l_NoCaspian.nc4"

# 3D monthly climatologies of potential temperature and salinity (can be those you used for the NEMO run):
export F_T_CLIM_3D_12="${CONF_INI_DIR}/thetao_1degx1deg-ORCA1.L75_WOA2009_monthly.nc4"
export F_S_CLIM_3D_12="${CONF_INI_DIR}/so_1degx1deg-ORCA1.L75_WOA2009_monthly.nc4"
export SST_CLIM_12="${STORE_DIR}/ORCA1.L75/ORCA1.L75-I/tos_180x360-ORCA1_Reynolds_monthly_mean1982-2005.nc"
export NN_T_CLIM="thetao"
export NN_S_CLIM="so"
export NN_SST_CLIM="tos"

export ICE_CLIM_12="${STORE_DIR}/ORCA1.L75/ORCA1.L75-I/ice_cover_180x360-ORCA1_Hurrell_monthly_mean1980-1999.nc4"
export NN_ICEF_CLIM="ice_cover"


# A text file where the vertical hydraugraphical sections of interest are defined :
export TRANSPORT_SECTION_FILE="${BARAKUDA_ROOT}/data/transportiz_ORCA1.dat"

# For transport by sigma-class:
export DENSITY_SECTION_FILE="${BARAKUDA_ROOT}/data/dens_section_ORCA1.dat"

# Files with the list of rectangular domains to "analyze" more closely:
export FILE_DEF_BOXES="${BARAKUDA_ROOT}/data/def_boxes_convection_ORCA1.txt"
export FILE_DMV_BOXES="${BARAKUDA_ROOT}/data/def_boxes_convection_ORCA1.txt"

# In what format should figures be produced ('png' recommanded, but 'svg' supported!):
export FIG_FORM="png"

# About remote HOST to install HTML pages to:
export ihttp=1 ; # do we export on a remote http server (1) or keep on the local machine (0)
export RHOST=misu228.misu.su.se ; # remote host to send diagnostic page to///
export RUSER=laurent ; # username associated to remote host (for file export)
export RWWWD=/data/www/barakuda/CMIP5 ; # directory of the local or remote host to send the diagnostic page to



#########################
# Diags to be performed #
#########################

# Movies of SST and SSS compared to OBS:
export i_do_movi=0

# Basic 3D and surface averages:
export i_do_mean=0

# IFS surface fluxes of heat and freshwater
export i_do_ifs_flx=0 ; # only relevant when ece_run=2...

# AMOC:
export i_do_amoc=0
export LMOCLAT="20-23 30-33 40-43 45-48 50-53" ; # List of latitude bands to look in for max of AMOC

# Transport of mass, heat and salt through specified sections (into TRANSPORT_SECTION_FILE):
export i_do_trsp=0  ; # transport of mass, heat and salt through specified sections
#              # i_do_trsp=2 => treat also different depths range!
z1_trsp=100  ; # first  depth: i_do_trsp must be set to 2
z2_trsp=1000 ; # second depth: i_do_trsp must be set to 2

# Meridional heat/salt transport (advective)
export i_do_mht=0

# Transport by sigma class
export i_do_sigt=0

# Sea-ice diags
export i_do_ice=0  ; # Sea-ice diags

# Budget on pre-defined (FILE_DEF_BOXES) rectangular domains:
export i_do_bb=0   ; # Budget and other stuffs on a given rectangular box!
#             # => needs file FILE_DEF_BOXES !!!
# => produces time-series f(t)  (mean of 2D fields)

# Vertical profiles on of box-averaged as a function of time...
export i_do_box_TS_z=0 ; # do sigma vert. profiles on given boxes... # 1 => no figures, 2 => figures
#                 # => needs file FILE_DEF_BOXES !!!
# => produces time-series f(t,z)

# Deep Mixed volume in prescribed boxes:
export i_do_dmv=0
export MLD_CRIT="1000,725,500"

# User-defined meridional or zonal cross sections (for temperature and salinity)
# => TS_SECTION_FILE must be defined!
export i_do_sect=1
export TS_SECTION_FILE="${BARAKUDA_ROOT}/data/TS_sections.dat"



# BETA / TESTING / NERDY (at your own risks...):
#
export i_do_ssx_box=0 ; # zoom on given boxes (+spatially-averaged values) for surface properties
#                     # boxes defined into barakuda_orca.py ...

# Some nerdy stuffs about the critical depth in prescribed boxes:
export i_do_zcrit=0

# Fresh-water transport associated to sea-ice transport
#  => must compile cdficeflux.x but depends on more recent CDFTOOLS module...
export i_do_icet=0 ; # treat sea-ice volume transport!
export TRANSPORT_ICE_SECTION_FILE="${BARAKUDA_ROOT}/data/transportiz_ORCA1_ARCTIC.dat"

export i_do_amo=0 ;  # buit a SST time serie usable to build Atlantic Multidecadal Oscilation index




# Place for potential specific host-related survival tricks:

#========================== Marenostrum @ BSC =========================================================
### Shouldn't be needed elsewhere than MareNostrum, where it's a hello to have CDO working...
## => Only if you specified ece_run=2 and i_do_ifs_flx
#export MOD_CDO="gcc/4.7.2 intel/13.0.1 openmpi/1.8.1 NETCDF/4.1.3 HDF5/1.8.10 UDUNITS/2.1.24 CDO/1.7.0"
#=======================================================================================================

