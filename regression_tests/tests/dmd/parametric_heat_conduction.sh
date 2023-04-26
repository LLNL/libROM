#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "rm -rf parameters.txt" 
    "$COMMAND ./parametric_heat_conduction -r 0.4  -offline -rdim 16" 
    "$COMMAND ./parametric_heat_conduction -r 0.45 -offline -rdim 16"
    "$COMMAND ./parametric_heat_conduction -r 0.55 -offline -rdim 16"
    "$COMMAND ./parametric_heat_conduction -r 0.6 -offline -rdim 16"
    "$COMMAND ./parametric_heat_conduction -r 0.5 -online -predict"
)
TYPE="DMD"
OFFSET=5
run_tests