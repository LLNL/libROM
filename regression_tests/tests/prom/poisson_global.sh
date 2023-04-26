#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "./poisson_global_rom -offline -f 1.0 -id 0" 
    "./poisson_global_rom -offline -f 1.1 -id 1" 
    "./poisson_global_rom -offline -f 1.2 -id 2" 
    "./poisson_global_rom -merge -ns 3" 
    "./poisson_global_rom -fom -f 1.15"
    "./poisson_global_rom -online -f 1.15"
)
TYPE="PROM"
OFFSET=5
run_tests
