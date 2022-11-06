#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "./dg_advection_global_rom -offline -ff 1.0 -id 0" 
    "./dg_advection_global_rom -offline -ff 1.1 -id 1" 
    "./dg_advection_global_rom -offline -ff 1.2 -id 2" 
    "./dg_advection_global_rom -merge -ns 3" 
    "./dg_advection_global_rom -online -ff 1.2"
)
TYPE="PROM"
run_tests
#cd ${EX_PROM_PATH_LOCAL}
#run_cmds

#cd ${EX_PROM_PATH_BASELINE} # Baseline(master) branch libROM
#run_cmds

#cd ${GITHUB_WORKSPACE}/build/tests

#if [[ ! -z test_offline ]]; then
   # ./basisComparator ${EX_PROM_PATH_LOCAL}/basis ${EX_PROM_PATH_BASELINE}/basis 1e-7 1
   # check_fail
#fi

#./solutionComparator ${EX_PROM_PATH_LOCAL}/dg_advection_global_rom-final.1.000000  ${EX_PROM_PATH_BASELINE}/dg_advection_global_rom-final.1.000000 "1.0e-5" "1"
#check_fail
#./solutionComparator ${EX_PROM_PATH_LOCAL}/dg_advection_global_rom-final.1.1.000000  ${EX_PROM_PATH_BASELINE}/dg_advection_global_rom-final.1.1.000000 "1.0e-5" "1"
#check_fail
#./solutionComparator ${EX_PROM_PATH_LOCAL}/dg_advection_global_rom-final.1.2.000000  ${EX_PROM_PATH_BASELINE}/dg_advection_global_rom-final.1.2.000000 "1.0e-5" "1"
#check_fail

#move_output_files










