#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "$COMMAND ./dg_euler -p 1 -rs 1 -rp 1 -o 5 -s 6 -tf 0.1" 
)
TYPE="DMD"
OFFSET=5
run_tests