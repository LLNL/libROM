#!/bin/bash
set -eo pipefail
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "mpirun -np 8 nonlinear_elasticity -s 2 -rs 1 -dt 0.01 -tf 5" 
)
TYPE="DMD"
cd ${EX_DMD_PATH_LOCAL}
run_cmds
cd ${EX_DMD_PATH_BASELINE} # Baseline(master) branch libROM
run_cmds

cd ${GITHUB_WORKSPACE}/build/tests

./solutionComparator ${EX_DMD_PATH_LOCAL}/elastic_energy.000000 ${EX_DMD_PATH_BASELINE}/elastic_energy.000000 "1.0e-5" "$1"

move_output_files









