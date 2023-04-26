#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "./nonlinear_elasticity_global_rom --offline -dt 0.01 -tf 5.0 -s 14 -vs 100 -sc 0.90 -id 0" 
    "./nonlinear_elasticity_global_rom --offline -dt 0.01 -tf 5.0 -s 14 -vs 100 -sc 1.10 -id 1" 
    "./nonlinear_elasticity_global_rom --merge -ns 2 -dt 0.01 -tf 5.0" 
    "./nonlinear_elasticity_global_rom -offline -dt 0.01 -tf 5.0 -s 14 -vs 100 -sc 1.00 -id 2"
    "./nonlinear_elasticity_global_rom --online -dt 0.01 -tf 5.0 -s 14 -vs 100 -hyp -rvdim 40 -rxdim 10 -hdim 71 -nsr 1170 -sc 1.00"
)
TYPE="PROM"
OFFSET=0
run_tests
