#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh
cd ${EX_PROM_PATH_LOCAL}
./dg_advection_global_rom -offline -ff 1.0 -id 0
./dg_advection_global_rom -offline -ff 1.1 -id 1
./dg_advection_global_rom -offline -ff 1.2 -id 2
./dg_advection_global_rom -merge -ns 3


cd ${EX_PROM_PATH_BASELINE} # Baseline(master) branch libROM
./dg_advection_global_rom -offline -ff 1.0 -id 0
./dg_advection_global_rom -offline -ff 1.1 -id 1
./dg_advection_global_rom -offline -ff 1.2 -id 2
./dg_advection_global_rom -merge -ns 3

cd ${GITHUB_WORKSPACE}/build/tests

./solutionComparator ${EX_PROM_PATH_LOCAL}/dg_advection_global_rom-final.1.000000  ${EX_PROM_PATH_BASELINE}/dg_advection_global_rom-final.1.000000 "1.0e-5" "1"
check_fail












