#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh
CMDS=( 
    "./mixed_nonlinear_diffusion -offline"
    "./mixed_nonlinear_diffusion -merge -ns 1"
    "./mixed_nonlinear_diffusion -online -rrdim 8 -rwdim 8 -nldim 20 -nsdim 20"
)
TYPE="PROM"
cd ${EX_PROM_PATH_LOCAL}
run_cmds

cd ${EX_PROM_PATH_BASELINE} # Baseline(master) branch libROM
run_cmds

cd ${GITHUB_WORKSPACE}/build/tests

./basisComparator ${EX_PROM_PATH_LOCAL}/basisR.000000 ${EX_PROM_PATH_BASELINE}/basisR.000000 1e-7 1
check_fail

./basisComparator ${EX_PROM_PATH_LOCAL}/basisFR.000000 ${EX_PROM_PATH_BASELINE}/basisFR.000000 1e-7 1
check_fail

./basisComparator ${EX_PROM_PATH_LOCAL}/basisW.000000 ${EX_PROM_PATH_BASELINE}/basisW.000000 1e-7 1
check_fail

./basisComparator ${EX_PROM_PATH_LOCAL}/basisS.000000 ${EX_PROM_PATH_BASELINE}/basisS.000000 1e-7 1
check_fail

./solutionComparator ${EX_PROM_PATH_LOCAL}/nldiff-final0.000000 ${EX_PROM_PATH_BASELINE}/nldiff-final0.000000 "1.0e-5" "1" 
check_fail

./solutionComparator ${EX_PROM_PATH_LOCAL}/nldiff-rom-final0.000000 ${EX_PROM_PATH_BASELINE}/nldiff-rom-final0.000000 "1.0e-5" "1" 
check_fail

move_output_files











