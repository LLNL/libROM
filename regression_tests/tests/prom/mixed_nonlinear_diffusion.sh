#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh
CMDS=( 
    "./mixed_nonlinear_diffusion -offline -tf 0.01"
    "./mixed_nonlinear_diffusion -merge -ns 1 -tf 0.01"
    "./mixed_nonlinear_diffusion -online -tf 0.01"
)
TYPE="PROM"
OFFSET=5
run_tests