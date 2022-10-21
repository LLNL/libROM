#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

BASELINE_DIR=$GITHUB_WORKSPACE/dependencies

CMDS=( 
    "rm -rf parameters.txt"
    "rm -rf de_parametric_heat_conduction_greedy_*"
    "$COMMAND de_parametric_heat_conduction_greedy -build_database -rdim 16 -greedy-param-size 5 -greedysubsize 2 -greedyconvsize 3 -greedyreldifftol 0.01" 
    "$COMMAND de_parametric_heat_conduction_greedy -r 0.2 -cx 0.2 -cy 0.2"
    "$COMMAND de_parametric_heat_conduction_greedy -r 0.2 -cx 0.2 -cy 0.2 -de -de_f 0.9 -de_cr 0.9 -de_ps 50 -de_min_iter 10 -de_max_iter 100 -de_ct 0.001 "
)
TYPE="DMD"
cd ${EX_DMD_PATH_LOCAL}
run_cmds

cd ${EX_DMD_PATH_BASELINE} 
run_cmds

cd ${GITHUB_WORKSPACE}/build/tests


./solutionComparator ${EX_DMD_PATH_LOCAL}/de_parametric_heat_conduction_greedy_0.105262_0.100000_0.141383_0.274614-final.000000  ${EX_DMD_PATH_BASELINE}/de_parametric_heat_conduction_greedy_0.105262_0.100000_0.141383_0.274614-final.000000  "1.0e-5" "8"
check_fail

./solutionComparator ${EX_DMD_PATH_LOCAL}/de_parametric_heat_conduction_greedy_0.106331_0.100000_0.107012_0.121436-final.000000  ${EX_DMD_PATH_BASELINE}/de_parametric_heat_conduction_greedy_0.106331_0.100000_0.107012_0.121436-final.000000  "1.0e-5" "8"
check_fail

./solutionComparator ${EX_DMD_PATH_LOCAL}/de_parametric_heat_conduction_greedy_0.138397_0.100000_0.128142_0.144331-final.000000  ${EX_DMD_PATH_BASELINE}/de_parametric_heat_conduction_greedy_0.138397_0.100000_0.128142_0.144331-final.000000 "1.0e-5" "8"
check_fail

./solutionComparator ${EX_DMD_PATH_LOCAL}/de_parametric_heat_conduction_greedy_0.106331_0.100000_0.107012_0.121436-final.000000  ${EX_DMD_PATH_BASELINE}/de_parametric_heat_conduction_greedy_0.106331_0.100000_0.107012_0.121436-final.000000 "1.0e-5" "8"
check_fail

./solutionComparator ${EX_DMD_PATH_LOCAL}/de_parametric_heat_conduction_greedy_0.180308_0.100000_0.207471_0.217217-final.000000  ${EX_DMD_PATH_BASELINE}/de_parametric_heat_conduction_greedy_0.180308_0.100000_0.207471_0.217217-final.000000 "1.0e-5" "8"
check_fail

./solutionComparator ${EX_DMD_PATH_LOCAL}/de_parametric_heat_conduction_greedy_0.200000_0.010000_0.200000_0.200000-final.000000 ${EX_DMD_PATH_BASELINE}/de_parametric_heat_conduction_greedy_0.200000_0.010000_0.200000_0.200000-final.000000 "1.0e-5" "8"
check_fail

./solutionComparator ${EX_DMD_PATH_LOCAL}/de_parametric_heat_conduction_greedy_0.214872_0.100000_0.233410_0.209747-final.000000 ${EX_DMD_PATH_BASELINE}/de_parametric_heat_conduction_greedy_0.214872_0.100000_0.233410_0.209747-final.000000 "1.0e-5" "8"
check_fail
./solutionComparator ${EX_DMD_PATH_LOCAL}/de_parametric_heat_conduction_greedy_0.217203_0.100000_0.252469_0.115465-final.000000 ${EX_DMD_PATH_BASELINE}/de_parametric_heat_conduction_greedy_0.217203_0.100000_0.252469_0.115465-final.000000"1.0e-5" "8"
check_fail
./solutionComparator ${EX_DMD_PATH_LOCAL}/de_parametric_heat_conduction_greedy_0.233593_0.100000_0.188263_0.200455-final.000000 ${EX_DMD_PATH_BASELINE}/de_parametric_heat_conduction_greedy_0.233593_0.100000_0.188263_0.200455-final.000000 "1.0e-5" "8"
check_fail

move_output_files
