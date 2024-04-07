#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh
CMDS=( 
    "$COMMAND ./heat_conduction_incdmd --inc" 
)
TYPE="DMD"
OFFSET=5
run_tests
