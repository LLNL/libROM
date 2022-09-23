#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "rm -rf parameters.txt" 
    "$COMMAND -np 8 parametric_heat_conduction -r 0.4 -visit -offline -rdim 16" 
    "$COMMAND -np 8 parametric_heat_conduction -r 0.45 -visit -offline -rdim 16"
    "$COMMAND -np 8 parametric_heat_conduction -r 0.55 -visit -offline -rdim 16"
    "$COMMAND -np 8 parametric_heat_conduction -r 0.6 -visit -offline -rdim 16"
    "$COMMAND -np 8 parametric_heat_conduction -r 0.5 -visit -online -predict"
)
TYPE="DMD"

cd ${EX_DMD_PATH_LOCAL}
run_cmds


cd ${EX_DMD_PATH_BASELINE} # Baseline(master) branch libROM
run_cmds

cd ${GITHUB_WORKSPACE}/build/tests

#./solutionComparator ${EX_DMD_PATH_LOCAL}/heat_conduction-final.000000 ${EX_DMD_PATH_BASELINE}/heat_conduction-final.000000 "1.0e-5" "$1"
check_fail

move_output_files