#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

BASELINE_DIR=$GITHUB_WORKSPACE/dependencies

cd ${GITHUB_WORKSPACE}/build/examples/dmd
mpirun -np 8 dg_advection -p 0 -dt 0.01 -tf 2

cd ${BASELINE_DIR}/libROM/build/examples/dmd # Baseline(master) branch libROM
 mpirun -np 8 dg_advection -p 0 -dt 0.01 -tf 2

cd ${GITHUB_WORKSPACE}/build/tests

# ./basisComparator ${GITHUB_WORKSPACE}/build/examples/dmd/heat_conduction-final ${BASELINE_DIR}/libROM/build/examples/dmd/heat_conduction-final 1e-7 1
# ./fileComparator ${GITHUB_WORKSPACE}/build/examples/dmd/dg_advection-final.000000 ${BASELINE_DIR}/libROM/build/examples/dmd/dg_advection-final.000000 "1.0e-5"
# check_fail
./solutionComparator ${GITHUB_WORKSPACE}/build/examples/dmd/dg_advection-final.000000 ${BASELINE_DIR}/libROM/build/examples/dmd/dg_advection-final.000000 "1.0e-5" "$1"
check_fail













