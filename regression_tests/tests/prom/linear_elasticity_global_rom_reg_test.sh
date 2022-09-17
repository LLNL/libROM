#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh
cd ${EX_PROM_PATH_LOCAL}
./linear_elasticity_global_rom -offline -id 0 -nu 0.2
./linear_elasticity_global_rom -offline -id 1 -nu 0.4
./linear_elasticity_global_rom -merge -ns 2

cd ${EX_PROM_PATH_BASELINE} # Baseline(master) branch libROM
./linear_elasticity_global_rom -offline -id 0 -nu 0.2
./linear_elasticity_global_rom -offline -id 1 -nu 0.4
./linear_elasticity_global_rom -merge -ns 2

cd ${GITHUB_WORKSPACE}/build/tests

./basisComparator ${EX_PROM_PATH_LOCAL}/basis ${EX_PROM_PATH_BASELINE}/basis "1e-7" "$1" 2>&1
check_fail
./solutionComparator ${EX_PROM_PATH_LOCAL}/Example_linear_elastic_000000/solution.000000  ${EX_PROM_PATH_BASELINE}/Example_linear_elastic_000000/solution.000000 "1.0e-5" "1"
check_fail

