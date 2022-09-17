#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

cd ${EX_PROM_PATH_LOCAL}
 
./mixed_nonlinear_diffusion -offline

cd ${EX_PROM_PATH_BASELINE} # Baseline(master) branch libROM
./mixed_nonlinear_diffusion -offline

cd ${GITHUB_WORKSPACE}/build/tests


./solutionComparator ${EX_PROM_PATH_LOCAL}/Sol0 ${EX_PROM_PATH_BASELINE}/Sol0 "1.0e-5" "1" 
check_fail












