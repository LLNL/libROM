#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh
cd ${GITHUB_WORKSPACE}/build/examples/prom
./dg_advection_global_rom -offline -ff 1.0 -id 0
./dg_advection_global_rom -offline -ff 1.1 -id 1
./dg_advection_global_rom -offline -ff 1.2 -id 2
./dg_advection_global_rom -merge -ns 3


cd ${BASELINE_DIR}/libROM/build/examples/prom # Baseline(master) branch libROM
./dg_advection_global_rom -offline -ff 1.0 -id 0
./dg_advection_global_rom -offline -ff 1.1 -id 1
./dg_advection_global_rom -offline -ff 1.2 -id 2
./dg_advection_global_rom -merge -ns 3

cd ${GITHUB_WORKSPACE}/build/tests

./basisComparator ${GITHUB_WORKSPACE}/build/examples/prom/basis ${BASELINE_DIR}/libROM/build/examples/prom/basis "1e-7 1" "$1"
check_fail
./solutionComparator ${GITHUB_WORKSPACE}/build/examples/prom/dg_advection_global_rom-final.1.000000  ${BASELINE_DIR}/libROM/build/examples/prom/dg_advection_global_rom-final.1.000000 "1.0e-5" "1"
check_fail












