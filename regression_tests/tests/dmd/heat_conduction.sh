
#!/bin/bash

set_pass() {
	echo "Heat_Conduction_Regression_Test: PASS"		
}


set_fail(){
 
	"Heat_Conduction_Regression_Test: FAIL"	
	
}

echo "Running Heat_Conduction_Regression_Test"

echo "GITHUB_WORKSPACE=$GITHUB_WORKSPACE"
BASELINE_DIR=$GITHUB_WORKSPACE/dependencies

cd ${GITHUB_WORKSPACE}/build/examples/dmd
mpirun -np 8 heat_conduction -s 1 -a 0.0 -k 1.0 
mpirun -np 8 heat_conduction -s 3 -a 0.5 -k 0.5 -o 4 -tf 0.7 -vs 1


cd ${BASELINE_DIR}/libROM/build/examples/dmd # Baseline(master) branch libROM
mpirun -np 8 heat_conduction -s 1 -a 0.0 -k 1.0 
mpirun -np 8 heat_conduction -s 3 -a 0.5 -k 0.5 -o 4 -tf 0.7 -vs 1

#cd ${GITHUB_WORKSPACE}/build/tests

#./basisComparator ${GITHUB_WORKSPACE}/build/examples/prom/basis ${BASELINE_DIR}/libROM/build/examples/prom/basis 1e-7 1
#./solutionComparator ${GITHUB_WORKSPACE}/build/examples/prom/sol.000000 ${BASELINE_DIR}/libROM/build/examples/prom/sol_f-0.01_fom.000000 "1.0e-5" "1" 

#echo "PIPESTATUS[0] = ${PIPESTATUS[0]}"
if [[ "${PIPESTATUS[0]}" -ne 0 ]];  # Capture and output the pipe status from MPI_Abort in some way
then
    set_fail

else
    set_pass

fi












