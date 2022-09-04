#!/bin/bash

tests_to_execute=("poisson_global_reg_test.sh" "dg_advection_global_rom_reg_test.sh" "linear_elasticity_global_rom_reg_test.sh")

echo "Setting up test suite"
echo "For detailed logs of the regression tests, please check regression_tests/results."

export GITHUB_WORKSPACE=/Users/pranav/Core/ROM_dev/libROM
DEPENDENCIES_DIR=${GITHUB_WORKSPACE}/dependencies
BUILD_DIR=${DEPENDENCIES_DIR}/build
MYDIR=$(pwd)
echo "My current dir = $MYDIR"
if [ ! -d $DEPENDENCIES_DIR ]; then
   echo "Creating $DEPENDENCIES_DIR"
   mkdir -p $DEPENDENCIES_DIR
   cd ${DEPENDENCIES_DIR}
   echo "Clone libROM master into dependencies"
   git clone https://github.com/LLNL/libROM.git
   cd libROM/scripts
   echo "Compile libROM master"
   ./compile.sh -m
   echo "Compile libROM master - done"
else
   echo "${DEPENDENCIES_DIR} already exists"
fi
cd $MYDIR
RESULTS_DIR=$MYDIR/results
if [ ! -d $RESULTS_DIR ]; then
    echo "Creating $RESULTS_DIR"
    mkdir -p $RESULTS_DIR;
else
	rm -rf $RESULTS_DIR/*
fi

testNum=${#tests_to_execute[@]}
echo "testNum = $testNum"
testNumPass=0
testNumFail=0

for test in ${tests_to_execute[@]}; do 
     ./$test
    if [[ "${PIPESTATUS[0]}" -ne 0 ]];  
        then
          testNumFail=$((testNumFail+1))

        else
          testNumPass=$((testNumPass+1))
    fi
done


echo "${testNumPass} passed, ${testNumFail} failed out of ${testNum} tests"

if [[ $testNumFail -ne 0 ]]; then
	exit 1
fi








