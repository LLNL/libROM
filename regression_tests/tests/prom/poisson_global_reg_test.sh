
#!/bin/bash

set_pass() {
	echo "Poisson_Global_Regression_Test: PASS"		
}


set_fail(){
 
	"Poisson_Global_Regression_Test: FAIL"	
	
}

echo "Running Poisson_Global_Regression_Test"

echo "GITHUB_WORKSPACE=$GITHUB_WORKSPACE"
BASELINE_DIR=$GITHUB_WORKSPACE/dependencies

cd ${GITHUB_WORKSPACE}/build/examples/prom
./poisson_global_rom -offline -f 1.0 -id 0
./poisson_global_rom -offline -f 1.1 -id 1
./poisson_global_rom -offline -f 1.2 -id 2
./poisson_global_rom -merge -ns 3
./poisson_global_rom -online -f 1.15


cd ${BASELINE_DIR}/libROM/build/examples/prom # Baseline(master) branch libROM
./poisson_global_rom -offline -f 1.0 -id 0
./poisson_global_rom -offline -f 1.1 -id 1
./poisson_global_rom -offline -f 1.2 -id 2
./poisson_global_rom -merge -ns 3
./poisson_global_rom -online -f 1.15

cd ${GITHUB_WORKSPACE}/build/tests

./basisComparator ${GITHUB_WORKSPACE}/build/examples/prom/basis ${BASELINE_DIR}/libROM/build/examples/prom/basis 1e-7 1
./solutionComparator ${GITHUB_WORKSPACE}/build/examples/prom/sol ${BASELINE_DIR}/libROM/build/examples/prom/sol_f-0.01_fom "1.0e-5" "1" 

#echo "PIPESTATUS[0] = ${PIPESTATUS[0]}"
if [[ "${PIPESTATUS[0]}" -ne 0 ]];  # Capture and output the pipe status from MPI_Abort in some way
then
    set_fail

else
    set_pass

fi












