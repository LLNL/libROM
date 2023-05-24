#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh
CMDS=( 
    "./wave_equation -o 4 -tf 5 -nwinsamp 25" 
)
TYPE="DMD"
OFFSET=5
run_tests
