#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

cd ${EX_PROM_PATH_LOCAL}
rm -rf frequencies.txt
./dg_advection_local_rom_matrix_interp -offline -ff 1.02 

cd ${EX_PROM_PATH_BASELINE} # Baseline(master) branch libROM
./dg_advection_local_rom_matrix_interp -offline -ff 1.02 

cd ${GITHUB_WORKSPACE}/build/tests

./solutionComparator ${EX_PROM_PATH_LOCAL}/Sol0 ${EX_PROM_PATH_BASELINE}/Sol0 "1.0e-5" "1" 
check_fail












