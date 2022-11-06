#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "./nonlinear_elasticity_global_rom --offline -dt 0.01 -tf 5.0 -s 14 -vs 100 -sc 0.90 -id 0" 
    "./nonlinear_elasticity_global_rom --offline -dt 0.01 -tf 5.0 -s 14 -vs 100 -sc 1.10 -id 1" 
    "./nonlinear_elasticity_global_rom --merge -ns 2 -dt 0.01 -tf 5.0" 
    "./nonlinear_elasticity_global_rom --online -dt 0.01 -tf 5.0 -s 14 -vs 100 -hyp -rvdim 40 -rxdim 10 -hdim 71 -nsr 1170 -sc 1.00"
)
TYPE="PROM"
run_tests
#cd ${EX_PROM_PATH_LOCAL}
#run_cmds

#cd ${EX_PROM_PATH_BASELINE} # Baseline(master) branch libROM
#run_cmds

#cd ${GITHUB_WORKSPACE}/build/tests

#./basisComparator ${EX_PROM_PATH_LOCAL}/basisH ${EX_PROM_PATH_BASELINE}/basisH 1e-7 1
#check_fail

#./basisComparator ${EX_PROM_PATH_LOCAL}/basisV ${EX_PROM_PATH_BASELINE}/basisV 1e-7 1
#check_fail

#./basisComparator ${EX_PROM_PATH_LOCAL}/basisX ${EX_PROM_PATH_BASELINE}/basisX 1e-7 1
#check_fail

#move_output_files










