
#!/bin/bash

RESULTS_DIR=$DIR/results
BASELINE_DIR=$GITHUB_WORKSPACE/dependencies

scriptName="dg_advection_global"
testName=""

#simulationLogFile="${RESULTS_DIR}/${scriptName}-${testName}.log"
#touch $simulationLogFile

set_pass() {
    echo "dg-advection-global-offline: PASS"		
}

set_fail(){
	echo "dg-advection-global-offline: FAIL"
}

cd ${GITHUB_WORKSPACE}/build/examples/prom
./dg_advection_global_rom -offline -ff 1.0 -id 0
./dg_advection_global_rom -offline -ff 1.1 -id 1
./dg_advection_global_rom -offline -ff 1.2 -id 2
./dg_advection_global_rom -merge -ns 3


cd ${BASELINE_DIR}/libROM/build/examples/prom # Baseline(master) branch libROM
./dg_advection_global_rom -offline -ff 1.0 -id 0
./dg_advection_global_rom -offline -ff 1.1 -id 1
./dg_advection_global_rom -offline -ff 1.2 -id 2
./dg_advection_global_rom -merge -ns 3

cd ${GITHUB_WORKSPACE}/build/tests

./basisComparator ${GITHUB_WORKSPACE}/build/examples/prom/basis ${BASELINE_DIR}/libROM/build/examples/prom/basis 1e-7 1

if [[ "${PIPESTATUS[0]}" -ne 0 ]];  # Capture and the pipe status from MPI_Abort to account for test failure
then
    set_fail

else
    set_pass

fi












