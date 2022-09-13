#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh
cd ${GITHUB_WORKSPACE}/build/examples/prom
./linear_elasticity_global_rom -offline -id 0 -nu 0.2
./linear_elasticity_global_rom -offline -id 1 -nu 0.4
./linear_elasticity_global_rom -merge -ns 2

cd ${BASELINE_DIR}/libROM/build/examples/prom # Baseline(master) branch libROM
./linear_elasticity_global_rom -offline -id 0 -nu 0.2
./linear_elasticity_global_rom -offline -id 1 -nu 0.4
./linear_elasticity_global_rom -merge -ns 2

cd ${GITHUB_WORKSPACE}/build/tests

./basisComparator ${GITHUB_WORKSPACE}/build/examples/prom/basis ${BASELINE_DIR}/libROM/build/examples/prom/basis "1e-7" "$1" 2>&1
check_fail
./solutionComparator ${GITHUB_WORKSPACE}/build/examples/prom/Example_linear_elastic_000000/solution.000000  ${BASELINE_DIR}/libROM/build/examples/prom/Example_linear_elastic_000000/solution.000000 "1.0e-5" "1"
check_fail

