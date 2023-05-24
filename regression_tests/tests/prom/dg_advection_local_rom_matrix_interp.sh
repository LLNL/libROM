#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "rm -rf frequencies.txt" 
    "./dg_advection_local_rom_matrix_interp -offline -ff 1.0"
    "./dg_advection_local_rom_matrix_interp -interp_prep -ff 1.0 -rdim 40"
    "./dg_advection_local_rom_matrix_interp -offline -ff 1.1"
    "./dg_advection_local_rom_matrix_interp -interp_prep -ff 1.1 -rdim 40"
    "./dg_advection_local_rom_matrix_interp -offline -ff 1.2"
    "./dg_advection_local_rom_matrix_interp -interp_prep -ff 1.2 -rdim 40"
    "./dg_advection_local_rom_matrix_interp -fom -ff 1.15"
    "./dg_advection_local_rom_matrix_interp -online_interp -ff 1.15 -rdim 40"
)
TYPE="PROM"
OFFSET=0
run_tests
