#!/bin/bash
set -eo pipefail
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "$COMMAND ./nonlinear_elasticity -s 2 -rs 1 -dt 0.01 -tf 5" 
)
TYPE="DMD"
OFFSET=5
run_tests
