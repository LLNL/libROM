
#!/bin/bash

set_pass() {
	echo "Poisson_Global_Regression_Test: PASS"		
}


set_fail(){
 
	"Poisson_Global_Regression_Test: FAIL"	
	
}

echo "Running Poisson_Global_Regression_Test"

echo "GITHUB_WORKSPACE=$GITHUB_WORKSPACE"

cd ${GITHUB_WORKSPACE}/build/examples/prom
./poisson_global_rom -offline -f 1.0 -id 0
./poisson_global_rom -offline -f 1.1 -id 1
./poisson_global_rom -offline -f 1.2 -id 2
./poisson_global_rom -merge -ns 3


cd ${GITHUB_WORKSPACE}/dependencies/libROM/build/examples/prom # Baseline(master) branch libROM
./poisson_global_rom -offline -f 1.0 -id 0
./poisson_global_rom -offline -f 1.1 -id 1
./poisson_global_rom -offline -f 1.2 -id 2
./poisson_global_rom -merge -ns 3

cd ${GITHUB_WORKSPACE}/build/tests

./basisComparator ${GITHUB_WORKSPACE}/build/examples/prom/basis ${GITHUB_WORKSPACE}/dependencies/libROM/build/examples/prom/basis 1e-7 1

echo "PIPESTATUS[0] = ${PIPESTATUS[0]}"
if [[ "${PIPESTATUS[0]}" -ne 0 ]];  # Capture and output the pipe status from MPI_Abort in some way
then
    set_fail

else
    set_pass

fi












