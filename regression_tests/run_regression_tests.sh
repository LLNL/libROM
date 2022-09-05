#!/bin/bash

# tests_to_execute=("poisson_global_reg_test.sh" "dg_advection_global_rom_reg_test.sh" "linear_elasticity_global_rom_reg_test.sh")

echo "Setting up test suite"
echo "For detailed logs of the regression tests, please check regression_tests/results."

export GITHUB_WORKSPACE=/Users/pranav/Core/ROM_dev/libROM
BASELINE_DIR=${GITHUB_WORKSPACE}/baseline
TESTS_DIR=${GITHUB_WORKSPACE}/regression_tests/tests
BUILD_DIR=${BASELINE_DIR}/build
MYDIR=$(pwd)
cd $TESTS_DIR
tests_to_execute=(*)
scriptName="Unknown"
echo "tests_to_execute = ${tests_to_execute}"
#echo "My current dir = $MYDIR"
if [ ! -d $BASELINE_DIR ]; then # Clone master branch to baseline directory
   echo "Creating $BASELINE_DIR"
   mkdir -p $BASELINE_DIR
   cd ${BASELINE_DIR}
   echo "Clone libROM master into dependencies"
   git clone https://github.com/LLNL/libROM.git
   cd libROM/scripts
   echo "Compile libROM master"
   ./compile.sh -m
   echo "Compile libROM master - done"
else
   echo "${BASELINE_DIR} already exists"
fi
cd $MYDIR
RESULTS_DIR=$MYDIR/results
if [ ! -d $RESULTS_DIR ]; then
    echo "Creating $RESULTS_DIR"
    mkdir -p $RESULTS_DIR
else
	rm -rf $RESULTS_DIR/*
fi

totalTests=${#tests_to_execute[@]}
echo "Number of tests = $totalTests"
testNum=0
testNumPass=0
testNumFail=0
rm -rf ${RESULTS_DIR}/* # Remove all log files from previous run

for test in ${tests_to_execute[@]}; do 
     scriptName=$(basename $test)
     echo "scriptName = $scriptName"
     simulationLogFile="${RESULTS_DIR}/${scriptName}.log"
     touch simulationLogFile
     testNum=$((testNum+1))
     ./tests/$test >> $simulationLogFile
    if [[ "${PIPESTATUS[0]}" -ne 0 ]];  
        then
          testNumFail=$((testNumFail+1))
          echo "$testNum. $test: FAIL"   
        else
          testNumPass=$((testNumPass+1))
          echo "$testNum. $test: PASS" 
          echo 
    fi
done


echo "${testNumPass} passed, ${testNumFail} failed out of ${totalTests} tests"
if [[ $testNumFail -ne 0 ]]; then
	exit 1
fi








