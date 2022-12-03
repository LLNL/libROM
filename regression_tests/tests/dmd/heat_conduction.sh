#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh
CMDS=( 
    "$COMMAND ./heat_conduction -s 1 -a 0.0 -k 1.0" 
)
TYPE="DMD"
OFFSET=5
run_tests
