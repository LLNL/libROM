#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh
CMDS=( 
    "./mixed_nonlinear_diffusion -offline -tf 0.01"
    "./mixed_nonlinear_diffusion -merge -ns 1 -tf 0.01"
    "./mixed_nonlinear_diffusion -online -rrdim 8 -rwdim 8 -nldim 20 -nsdim 20 -tf 0.01"
)
TYPE="PROM"
cd ${EX_PROM_PATH_LOCAL}
run_cmds

cd ${EX_PROM_PATH_BASELINE} # Baseline(master) branch libROM
run_cmds

cd ${GITHUB_WORKSPACE}/build/tests

./basisComparator ${EX_PROM_PATH_LOCAL}/basisR ${EX_PROM_PATH_BASELINE}/basisR 1e-7 1
check_fail

./basisComparator ${EX_PROM_PATH_LOCAL}/basisFR ${EX_PROM_PATH_BASELINE}/basisFR 1e-7 1
check_fail

./basisComparator ${EX_PROM_PATH_LOCAL}/basisW ${EX_PROM_PATH_BASELINE}/basisW 1e-7 1
check_fail

./basisComparator ${EX_PROM_PATH_LOCAL}/basisS ${EX_PROM_PATH_BASELINE}/basisS 1e-7 1
check_fail

./solutionComparator ${EX_PROM_PATH_LOCAL}/nldiff-final0.000000 ${EX_PROM_PATH_BASELINE}/nldiff-final0.000000 "1.0e-5" "1" 
check_fail

./solutionComparator ${EX_PROM_PATH_LOCAL}/nldiff-rom-final0.000000 ${EX_PROM_PATH_BASELINE}/nldiff-rom-final0.000000 "1.0e-5" "1" 
check_fail

move_output_files











