#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

cd  ${EX_PROM_PATH_LOCAL}
# build_database phase: 
./poisson_local_rom_greedy -build_database -greedy-param-min 1.0 -greedy-param-max 1.2 -greedy-param-size 5 -greedysubsize 2 -greedyconvsize 3 -greedyrelerrortol 0.01
# use_database phase:  (create a new solution to compare with) 
./poisson_local_rom_greedy -fom -f 1.15 
# use_database phase: (use the database to compute at f 1.15 while comparing to the true offline solution at f 1.15)   
./poisson_local_rom_greedy -use_database -online -f 1.15 

cd ${EX_PROM_PATH_BASELINE} # Baseline(master) branch libROM
./poisson_local_rom_greedy -build_database -greedy-param-min 1.0 -greedy-param-max 1.2 -greedy-param-size 5 -greedysubsize 2 -greedyconvsize 3 -greedyrelerrortol 0.01
./poisson_local_rom_greedy -fom -f 1.15 
./poisson_local_rom_greedy -use_database -online -f 1.15 

cd ${GITHUB_WORKSPACE}/build/tests

./solutionComparator ${EX_PROM_PATH_LOCAL}/Sol0 ${EX_PROM_PATH_BASELINE}/Sol0 "1.0e-5" "1" 
check_fail












