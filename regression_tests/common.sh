#!/bin/bash
# Common code for all the dmd and prom regression tests
set -eo pipefail
export RESULTS_DIR=$DIR/results
export BASELINE_DIR=$GITHUB_WORKSPACE/dependencies
export EX_DIR_LOCAL=${GITHUB_WORKSPACE}/build/examples
export EX_DMD_PATH_LOCAL=${EX_DIR_LOCAL}/dmd
export EX_PROM_PATH_LOCAL=${EX_DIR_LOCAL}/prom
export EX_DIR_BASELINE=${BASELINE_DIR}/libROM/build/examples
export EX_DMD_PATH_BASELINE=${EX_DIR_BASELINE}/dmd
export EX_PROM_PATH_BASELINE=${EX_DIR_BASELINE}/prom
export TYPE
trap "move_output_files_after_error" ERR
echo "GITHUB_WORKSPACE=$GITHUB_WORKSPACE"
echo "Running $0 with $1 processors"
SCRIPT_NAME=$(basename "$0" ".sh")
NP=$(($1))
# Check machine
echo "OS: $MACHINE"
if [[ $MACHINE = "Linux" ]]; then
	COMMAND="srun -p pbatch -n 8"
elif [[ $MACHINE = "Darwin" ]]; then
	COMMAND="mpirun -np 8"
else
    echo "Bad OS: $MACHINE"
    exit 1
fi

move_output_files_after_error() {
    echo "Moving output files after a trap"
    move_output_files
}

move_output_files() {
    EX_PATHS=("${EX_DMD_PATH_LOCAL}" "${EX_DMD_PATH_BASELINE}" "${EX_PROM_PATH_LOCAL}" "${EX_PROM_PATH_BASELINE}")
    echo "TYPE = $TYPE"
    echo "SCRIPT_NAME = $SCRIPT_NAME"
    for path in "${EX_PATHS[@]}"; do
        if [[ "$TYPE" == "DMD" && ($path ==  "${EX_PROM_PATH_LOCAL}" || $path == "${EX_PROM_PATH_BASELINE}") ]]; then
            continue
        elif [[ "$TYPE" = "PROM" && ($path ==  "${EX_DMD_PATH_LOCAL}" || $path == "${EX_DMD_PATH_BASELINE}") ]]; then
            continue
        fi
        cd "${path}"
        echo "Removing ${SCRIPT_NAME}_out in ${path}"
        rm -rf "${SCRIPT_NAME}_out"
        mkdir "${SCRIPT_NAME}_out"
        IFS=$'\n'
        MOVABLE_FILES=$(find . -maxdepth 1 ! -perm 755)
        unset $IFS
        
        for file in "${MOVABLE_FILES[@]}"; do
            echo "file = $file"
            re='_out'
            if [[ "${file}" =~ _out ]]; then
               echo "Skipping $file"
               continue
            fi
            echo "Moving $file to out directory"
            mv $file ${SCRIPT_NAME}_out
        done
    done
}

run_cmds() {
    for cmd in "${CMDS[@]}"; do
        eval "$cmd"
    done
}

#run_test() {

#}

set_pass() {
	echo "$SCRIPT_NAME: PASS"		
}

set_fail(){
	echo "$SCRIPT_NAME: FAIL"	
    exit 1
}

check_fail(){
    if [[ "$?" -ne 0 || "${PIPESTATUS[0]}" -ne 0 ]];  # Capture the pipe status from MPI_Abort 
    then
        set_fail
    else
        set_pass
    fi
}
