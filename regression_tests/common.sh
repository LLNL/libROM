#!/bin/bash
# Common code for all the dmd and prom regression tests
set -eo pipefail
RESULTS_DIR=$DIR/results
BASELINE_DIR=$GITHUB_WORKSPACE/dependencies
echo "GITHUB_WORKSPACE=$GITHUB_WORKSPACE"
echo "Running $0 with $1 processors"
SCRIPT_NAME=$(basename $0 ".sh")
NP=$(($1))
set_pass() {
	echo "$SCRIPT_NAME: PASS"		
}

set_fail(){
	echo "$SCRIPT_NAME: FAIL"	
    exit 1
}

check_fail(){
    if [[ "$?" -ne 0 || "${PIPESTATUS[0]}" -ne 0 ]];  # Capture and output the pipe status from MPI_Abort in some way
    then
        set_fail
    else
        set_pass
    fi
}
