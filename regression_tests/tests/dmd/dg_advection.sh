
#!/bin/bash

set_pass() {
	echo "dg_advection_Regression_Test: PASS"		
}


set_fail(){
	"dg_advection_Regression_Test: FAIL"	
}

echo "Running dg_advection regression test"

echo "GITHUB_WORKSPACE=$GITHUB_WORKSPACE"
BASELINE_DIR=$GITHUB_WORKSPACE/dependencies

cd ${GITHUB_WORKSPACE}/build/examples/dmd
mpirun -np 8 dg_advection -p 0 -dt 0.01 -tf 2
#mpirun -np 8 heat_conduction -s 3 -a 0.5 -k 0.5 -o 4 -tf 0.7 -vs 1


cd ${BASELINE_DIR}/libROM/build/examples/dmd # Baseline(master) branch libROM
 mpirun -np 8 dg_advection -p 0 -dt 0.01 -tf 2
#mpirun -np 8 heat_conduction -s 3 -a 0.5 -k 0.5 -o 4 -tf 0.7 -vs 1

cd ${GITHUB_WORKSPACE}/build/tests

echo "Running solution comparator"
#./basisComparator ${GITHUB_WORKSPACE}/build/examples/dmd/heat_conduction-final ${BASELINE_DIR}/libROM/build/examples/dmd/heat_conduction-final 1e-7 1
./solutionComparator ${GITHUB_WORKSPACE}/build/examples/dmd/dg_advection-final.000000 ${BASELINE_DIR}/libROM/build/examples/dmd/dg_advection-final.000000 "1.0e-5" "1" 

#echo "PIPESTATUS[0] = ${PIPESTATUS[0]}"
if [[ "${PIPESTATUS[0]}" -ne 0 ]];  # Capture and output the pipe status from MPI_Abort in some way
then
    set_fail

else
    set_pass

fi












