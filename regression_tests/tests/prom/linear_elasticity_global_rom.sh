#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "./linear_elasticity_global_rom -offline -id 0 -nu 0.2" 
    "./linear_elasticity_global_rom -offline -id 1 -nu 0.4" 
    "./linear_elasticity_global_rom -merge -ns 2"
    "./linear_elasticity_global_rom -online -id 3 -nu 0.3"
)
TYPE="PROM"
cd ${EX_PROM_PATH_LOCAL}
run_cmds

cd ${EX_PROM_PATH_BASELINE} # Baseline(master) branch libROM
run_cmds

cd ${GITHUB_WORKSPACE}/build/tests

if [[ ! -z test_offline ]]; then
    ./basisComparator ${EX_PROM_PATH_LOCAL}/basis ${EX_PROM_PATH_BASELINE}/basis "1e-7" "$1" 2>&1
    check_fail
fi
./solutionComparator ${EX_PROM_PATH_LOCAL}/sol_f-0.01_rom.000000  ${EX_PROM_PATH_BASELINE}/sol_f-0.01_rom.000000 "1.0e-5" "1"
check_fail

move_output_files
