#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "./dg_advection_global_rom -offline -ff 1.0 -id 0" 
    "./dg_advection_global_rom -offline -ff 1.1 -id 1" 
    "./dg_advection_global_rom -offline -ff 1.2 -id 2" 
    "./dg_advection_global_rom -merge -ns 3" 
    "./dg_advection_global_rom -fom -ff 1.15"
    "./dg_advection_global_rom -online -ff 1.15"
)
TYPE="PROM"
OFFSET=0
run_tests
