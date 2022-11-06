#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "$COMMAND dg_advection -p 0 -dt 0.01 -tf 2" 
)
TYPE="DMD"
run_tests
#cd ${EX_DMD_PATH_LOCAL}
#run_cmds

#cd ${EX_DMD_PATH_BASELINE} 
#run_cmds

#cd ${GITHUB_WORKSPACE}/build/tests

#./solutionComparator ${EX_DMD_PATH_LOCAL}/dg_advection-final.000000 ${EX_DMD_PATH_BASELINE}/dg_advection-final.000000 "1.0e-5" "8"
#check_fail

#move_output_files









