#!/bin/bash

tests_to_execute = ("poisson_global_reg_test.sh" "dg_advection_global_regression_test.sh" "linear_elasticity_global_rom_reg_test.sh")

echo "Setting up test suite"
echo "For detailed logs of the regression tests, please check regression_tests/results."

RESULTS_DIR=$DIR/results

if [ ! -d $RESULTS_DIR ]; then
    mkdir -p $RESULTS_DIR;
else
	rm -rf $RESULTS_DIR/*
fi

testNum=${#tests_to_execute[@]}
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








