#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "rm -rf frequencies.txt" 
    "./dg_advection_local_rom_matrix_interp --mesh ../data/periodic-square.mesh -offline -rs 4 -ff 1.02"
    "./dg_advection_local_rom_matrix_interp --mesh ../data/periodic-square.mesh-interp_prep -rs 4 -ff 1.02 -rdim 40"
    "./dg_advection_local_rom_matrix_interp --mesh ../data/periodic-square.mesh -offline -rs 4 -ff 1.08"
    "./dg_advection_local_rom_matrix_interp --mesh ../data/periodic-square.mesh -interp_prep -rs 4 -ff 1.08 -rdim 40"
    "./dg_advection_local_rom_matrix_interp --mesh ../data/periodic-square.mesh -fom -rs 4 -ff 1.05 -visit"
    "./dg_advection_local_rom_matrix_interp --mesh ../data/periodic-square.mesh -online_interp -rs 4 -ff 1.05 -rdim 40"
)
TYPE="PROM"
OFFSET=0
run_tests
