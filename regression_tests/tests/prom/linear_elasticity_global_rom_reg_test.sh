
#!/bin/bash

set_pass() {
	echo "linear_elasticity_global_rom_regression_test: PASS"	
    # echo "linear_elasticity_global_rom_regression_test: FAIL" >> $simulationLogFile
    #exit 1
}


set_fail(){
 
	echo "linear_elasticity_global_rom_regression_test: FAIL"
    # echo "linear_elasticity_global_rom_regression_test: FAIL" >> $simulationLogFile
   # exit 1

}

BASELINE_DIR=$GITHUB_WORKSPACE/dependencies

cd ${GITHUB_WORKSPACE}/build/examples/prom
./linear_elasticity_global_rom -offline -id 0 -nu 0.2
./linear_elasticity_global_rom -offline -id 1 -nu 0.4
./linear_elasticity_global_rom -merge -ns 2

cd ${BASELINE_DIR}/libROM/build/examples/prom # Baseline(master) branch libROM
./linear_elasticity_global_rom -offline -id 0 -nu 0.2
./linear_elasticity_global_rom -offline -id 1 -nu 0.4
./linear_elasticity_global_rom -merge -ns 2

cd ${GITHUB_WORKSPACE}/build/tests

./basisComparator ${GITHUB_WORKSPACE}/build/examples/prom/basis ${BASELINE_DIR}/libROM/build/examples/prom/basis 1e-7 1

if [[ "${PIPESTATUS[0]}" -ne 0 ]];  # Capture and output the pipe status from MPI_Abort
then
    set_fail

else
    set_pass

fi

