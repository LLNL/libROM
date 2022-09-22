#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh
CMDS=( 
    "./mixed_nonlinear_diffusion -offline"
)
TYPE="PROM"
cd ${EX_PROM_PATH_LOCAL}
run_cmds

cd ${EX_PROM_PATH_BASELINE} # Baseline(master) branch libROM
run_cmds

cd ${GITHUB_WORKSPACE}/build/tests


./solutionComparator ${EX_PROM_PATH_LOCAL}/nldiff-final0.000000 ${EX_PROM_PATH_BASELINE}/nldiff-final0.000000 "1.0e-5" "1" 
check_fail

move_output_files











