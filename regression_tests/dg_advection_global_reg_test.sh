
#!/bin/bash

RESULTS_DIR=$DIR/results

scriptName="dg_advection_global"
testName=""

simulationLogFile="${RESULTS_DIR}/${scriptName}-${testName}.log"
touch $simulationLogFile

set_pass() {
    echo "dg-advection-global-offline: PASS"		
    echo "dg-advection-global-offline: PASS" >> $simulationLogFile
}

set_fail(){
	echo "dg-advection-global-offline: FAIL"
    echo "dg-advection-global-offline: FAIL" >> $simulationLogFile
    exit 1
}

cd ${GITHUB_WORKSPACE}/build/examples/prom
./dg_advection_global_rom -offline -ff 1.0 -id 0
./dg_advection_global_rom -offline -ff 1.1 -id 1
./dg_advection_global_rom -offline -ff 1.2 -id 2
./dg_advection_global_rom -merge -ns 3


cd ${GITHUB_WORKSPACE}/dependencies/build/examples/prom # Baseline(master) branch libROM
./dg_advection_global_rom -offline -ff 1.0 -id 0
./dg_advection_global_rom -offline -ff 1.1 -id 1
./dg_advection_global_rom -offline -ff 1.2 -id 2
./dg_advection_global_rom -merge -ns 3

cd ${GITHUB_WORKSPACE}/build/tests/regression_tests

./basisComparator ${GITHUB_WORKSPACE}/build/examples/prom/basis.000000 ${GITHUB_WORKSPACE}/dependencies/build/examples/prom/basis.000000 1e-7 1


if [[ "${PIPESTATUS[0]}" -ne 0 ]];  # Capture and the pipe status from MPI_Abort to account for test failure
then
    set_fail

else
    set_pass

fi












