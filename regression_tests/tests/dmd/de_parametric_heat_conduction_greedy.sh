#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

BASELINE_DIR=$GITHUB_WORKSPACE/dependencies

CMDS=( 
    "rm -rf parameters.txt"
    "rm -rf de_parametric_heat_conduction_greedy_*"
    "$COMMAND ./de_parametric_heat_conduction_greedy -build_database -rdim 16 -greedy-param-size 5 -greedysubsize 2 -greedyconvsize 3 -greedyreldifftol 0.01" 
    "$COMMAND ./de_parametric_heat_conduction_greedy -r 0.2 -cx 0.2 -cy 0.2"
    "$COMMAND ./de_parametric_heat_conduction_greedy -r 0.2 -cx 0.2 -cy 0.2 -de -de_f 0.9 -de_cr 0.9 -de_ps 50 -de_min_iter 10 -de_max_iter 100 -de_ct 0.001 "
)
OFFSET=0
TYPE="DMD"
cd ${EX_DMD_PATH_LOCAL}
run_cmds
latest_file=$(ls -Art de_parametric_heat_conduction_greedy*final*|tail -n 1)
#echo "Latest file = $latest_file"

cd ${EX_DMD_PATH_BASELINE} 
run_cmds

cd ${GITHUB_WORKSPACE}/build/tests
./solutionComparator "${EX_DMD_PATH_BASELINE}/$latest_file" "${EX_DMD_PATH_LOCAL}/$latest_file" "1.0e-5" "1"
check_fail
move_output_files
