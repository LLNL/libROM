#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

cd ${EX_PROM_PATH_LOCAL}
./poisson_global_rom -offline -f 1.0 -id 0
./poisson_global_rom -offline -f 1.1 -id 1
./poisson_global_rom -offline -f 1.2 -id 2
./poisson_global_rom -merge -ns 3
./poisson_global_rom -online -f 1.15


cd ${EX_PROM_PATH_LOCAL} # Baseline(master) branch libROM
./poisson_global_rom -offline -f 1.0 -id 0
./poisson_global_rom -offline -f 1.1 -id 1
./poisson_global_rom -offline -f 1.2 -id 2
./poisson_global_rom -merge -ns 3
./poisson_global_rom -online -f 1.15

cd ${GITHUB_WORKSPACE}/build/tests

./basisComparator ${EX_PROM_PATH_LOCAL}/basis ${EX_PROM_PATH_BASELINE}/basis "1e-7" "1"
check_fail
./solutionComparator ${EX_PROM_PATH_LOCAL}/sol.000000 ${EX_PROM_PATH_BASELINE}/sol.000000 "1.0e-5" "1" 
check_fail












