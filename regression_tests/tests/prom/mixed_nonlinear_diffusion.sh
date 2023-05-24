#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh
CMDS=( 
    "./mixed_nonlinear_diffusion -offline -tf 0.01"
    "./mixed_nonlinear_diffusion -merge -ns 1 -tf 0.01"
    "./mixed_nonlinear_diffusion -online -tf 0.01 -rrdim 2 -rwdim 2 -nldim 3"
    "./mixed_nonlinear_diffusion -online -tf 0.01 -rrdim 2 -rwdim 2 -nldim 3 -sopt"
    "./mixed_nonlinear_diffusion -online -tf 0.01 -rrdim 2 -rwdim 2 -nldim 3 -ns 1 -eqp -maxnnls 4"
)
TYPE="PROM"
OFFSET=5
run_tests
