#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "./linear_elasticity_global_rom -offline -id 0 -nu 0.2" 
    "./linear_elasticity_global_rom -offline -id 1 -nu 0.4" 
    "./linear_elasticity_global_rom -merge -ns 2"
    "./linear_elasticity_global_rom -offline -id 2 -nu 0.3"
    "./linear_elasticity_global_rom -online -id 3 -nu 0.3"
)
TYPE="PROM"
OFFSET=0
run_tests
