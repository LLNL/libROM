#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

cd ${GITHUB_WORKSPACE}/build/examples/prom
./poisson_global_rom -offline -f 1.0 -id 0
./poisson_global_rom -offline -f 1.1 -id 1
./poisson_global_rom -offline -f 1.2 -id 2
./poisson_global_rom -merge -ns 3
./poisson_global_rom -online -f 1.15


cd ${BASELINE_DIR}/libROM/build/examples/prom # Baseline(master) branch libROM
./poisson_global_rom -offline -f 1.0 -id 0
./poisson_global_rom -offline -f 1.1 -id 1
./poisson_global_rom -offline -f 1.2 -id 2
./poisson_global_rom -merge -ns 3
./poisson_global_rom -online -f 1.15

cd ${GITHUB_WORKSPACE}/build/tests

./basisComparator ${GITHUB_WORKSPACE}/build/examples/prom/basis ${BASELINE_DIR}/libROM/build/examples/prom/basis "1e-7" "1"
check_fail
./solutionComparator ${GITHUB_WORKSPACE}/build/examples/prom/sol.000000 ${BASELINE_DIR}/libROM/build/examples/prom/sol.000000 "1.0e-5" "1" 
check_fail












