#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "$COMMAND ./dg_advection -p 0 -dt 0.01 -tf 2" 
)
TYPE="DMD"
OFFSET=5
run_tests
