#!/bin/bash

# tests_to_execute=("poisson_global_reg_test.sh" "dg_advection_global_rom_reg_test.sh" "linear_elasticity_global_rom_reg_test.sh")

echo "Setting up test suite"
echo "For detailed logs of the regression tests, please check regression_tests/results."
if [[ -z ${GITHUB_WORKSPACE} ]]; then
  echo "GITHUB_WORKSPACE is not set. Exiting ..."
  exit 1
fi
echo "GITHUB_WORKSPACE = ${GITHUB_WORKSPACE}"
export TMPDIR=/tmp
BASELINE_DIR=${GITHUB_WORKSPACE}/dependencies
TESTS_DIR=${GITHUB_WORKSPACE}/regression_tests/tests
BUILD_DIR=${BASELINE_DIR}/build
DIR=$(pwd)
scriptName="Unknown"
#echo "My current dir = $DIR"
if [ ! -d $BASELINE_DIR/libROM ]; then # Clone master branch to baseline directory
   echo "Creating $BASELINE_DIR"
   mkdir -p $BASELINE_DIR
   cd ${BASELINE_DIR}
   echo "Clone libROM master into baseline"
   git clone https://github.com/LLNL/libROM.git
   cd libROM/scripts
   echo "Compile libROM master"
   ./compile.sh -m
   echo "Compile libROM master - done"
else
   echo "${BASELINE_DIR}/libROM already exists"
fi
cd $DIR
RESULTS_DIR=$DIR/regression_tests/results
if [ ! -d $RESULTS_DIR ]; then
    echo "Creating $RESULTS_DIR"
    mkdir -p $RESULTS_DIR
else
  echo "Removing old files from $RESULTS_DIR"
	rm -rf $RESULTS_DIR/*
fi

# Get the number of processors
NUM_PROCESSORS=$(getconf _NPROCESSORS_ONLN)
re='^[0-9]+$'
if ! [[ $NUM_PROCESSORS =~ $re ]] ; then
   echo "Error: $NUM_PROCESSORS is not a number"
   exit 1
fi
echo "Number of processors = $NUM_PROCESSORS"
echo "Running regression tests"
totalTests=0
testNum=0
testNumPass=0
testNumFail=0
cd $TESTS_DIR
type_of_tests_to_execute=(*)
for type_of_test in ${type_of_tests_to_execute[@]}; do
  cd $type_of_test
  tests_to_execute=(*)
  for test in ${tests_to_execute[@]}; do 
      scriptName=$(basename $test ".sh")
      # Run a specific test by specifying the test (including the .sh suffix)
      if [[ -n "$1" && "$1" != "$test" ]]; then
         continue
      fi
          
      simulationLogFile="${RESULTS_DIR}/${scriptName}.log"
      touch $simulationLogFile
      testNum=$((testNum+1))
      ./$test "$NUM_PROCESSORS" >> $simulationLogFile 2>&1
      if [[ $? -ne 0 || "${PIPESTATUS[0]}" -ne 0 ]];  
          then
            testNumFail=$((testNumFail+1))
            echo "$testNum. $test: FAIL"   
          else
            testNumPass=$((testNumPass+1))
            echo "$testNum. $test: PASS" 
      fi
  done
  cd ..
done
totalTests=$testNum
echo "${testNumPass} passed, ${testNumFail} failed out of ${totalTests} tests"
if [[ $testNumFail -ne 0 ]]; then
	exit 1
fi








