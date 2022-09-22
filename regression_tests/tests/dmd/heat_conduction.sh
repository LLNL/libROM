#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "mpirun -np 8 heat_conduction -s 1 -a 0.0 -k 1.0" 
)
TYPE="DMD"
cd ${EX_DMD_PATH_LOCAL}
run_cmds

cd ${EX_DMD_PATH_BASELINE} # Baseline(master) branch libROM
run_cmds

cd ${GITHUB_WORKSPACE}/build/tests

./solutionComparator ${EX_DMD_PATH_LOCAL}/heat_conduction-final.000000 ${EX_DMD_PATH_BASELINE}/heat_conduction-final.000000 "1.0e-5" "$1"
check_fail

move_output_files










