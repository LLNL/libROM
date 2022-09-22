#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "rm -rf frequencies.txt" 
    "./dg_advection_local_rom_matrix_interp -offline -ff 1.02 "
)
TYPE="PROM"
cd ${EX_PROM_PATH_LOCAL}
run_cmds

cd ${EX_PROM_PATH_BASELINE} # Baseline(master) branch libROM
run_cmds

cd ${GITHUB_WORKSPACE}/build/tests

./basisComparator ${EX_PROM_PATH_LOCAL}/basis_1.020000 ${EX_PROM_PATH_BASELINE}/basis_1.020000 1e-7 1
check_fail

./solutionComparator ${EX_PROM_PATH_LOCAL}/dg_advection_local_rom_matrix_interp-final.1.02.000000 ${EX_PROM_PATH_BASELINE}/dg_advection_local_rom_matrix_interp-final.1.02.000000 "1.0e-5" "1" 
check_fail

move_output_files











