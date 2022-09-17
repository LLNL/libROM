#!/bin/bash
set -eo pipefail
source $GITHUB_WORKSPACE/regression_tests/common.sh

cd ${EX_DMD_PATH_LOCAL}
mpirun -np 8 nonlinear_elasticity -s 2 -rs 1 -dt 0.01 -tf 5

cd ${EX_DMD_PATH_BASELINE} # Baseline(master) branch libROM
mpirun -np 8 nonlinear_elasticity -s 2 -rs 1 -dt 0.01 -tf 5

cd ${GITHUB_WORKSPACE}/build/tests

#./basisComparator ${GITHUB_WORKSPACE}/build/examples/dmd/heat_conduction-final ${BASELINE_DIR}/libROM/build/examples/dmd/heat_conduction-final 1e-7 1
./solutionComparator ${EX_DMD_PATH_LOCAL}/elastic_energy.000000 ${EX_DMD_PATH_BASELINE}/elastic_energy.000000 "1.0e-5" "$1"











