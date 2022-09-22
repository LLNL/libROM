#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "./poisson_local_rom_greedy -build_database -greedy-param-min 1.0 -greedy-param-max 1.2 -greedy-param-size 5 -greedysubsize 2 -greedyconvsize 3 -greedyrelerrortol 0.01"  
    "./poisson_local_rom_greedy -fom -f 1.15" 
    "./poisson_local_rom_greedy -use_database -online -f 1.15" 
)
TYPE="PROM"
cd  ${EX_PROM_PATH_LOCAL}
run_cmds

cd ${EX_PROM_PATH_BASELINE} # Baseline(master) branch libROM
run_cmds

cd ${GITHUB_WORKSPACE}/build/tests

./solutionComparator ${EX_PROM_PATH_LOCAL}/Sol0 ${EX_PROM_PATH_BASELINE}/Sol0 "1.0e-5" "1" 
check_fail

move_output_files











