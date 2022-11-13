#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "rm -rf frequencies.txt" 
    "./dg_advection_local_rom_matrix_interp -offline"
    "./dg_advection_local_rom_matrix_interp -online"
)
TYPE="PROM"
run_tests











