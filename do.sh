#!/bin/bash

if [ "$1" = "e" ]; then
    sbatch -n 6 -t 1800 ./barakuda.sh -C ORCA025_L75_etienne -R a0ez -e
    exit
fi

sbatch -n 7 -t 36000 ./barakuda.sh -C ORCA025_L75_etienne -R a0ez 
