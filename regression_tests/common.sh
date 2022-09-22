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
trap "move_output_files" ERR
echo "GITHUB_WORKSPACE=$GITHUB_WORKSPACE"
echo "Running $0 with $1 processors"
SCRIPT_NAME=$(basename "$0" ".sh")
NP=$(($1))
DMD=false
PROM=false
if [[ -z "${PROGRAM_FILES}" ]]; then
    rm -rf EX_DMD_PATH_LOCAL/*/
    rm -rf EX_PROM_PATH_LOCAL/*/
    rm -rf EX_DMD_PATH_BASELINE/*/
    rm -rf EX_PROM_PATH_BASELINE/*/
    cd $EX_DMD_PATH_BASELINE
    DMD_PROGRAM_FILES=(*)
    echo "DMD_PROGRAM_FILES = $DMD_PROGRAM_FILES"
    cd $EX_PROM_PATH_BASELINE
    PROM_PROGRAM_FILES=(*)
    echo "PROM_PROGRAM_FIES = $PROM_PROGRAM_FILES"
    PROGRAM_FILES_WITH_PATH=("${DMD_PROGRAM_FILES[@]}" "${PROM_PROGRAM_FILES[@]}")
    echo "PROGRAM_FILES_WITH_PATH = $PROGRAM_FILES_WITH_PATH"
    PROGRAM_FILES=()
    for x in ${PROGRAM_FILES_WITH_PATH[@]}; do
        name=($(basename $x))
        PROGRAM_FILES+=("${name}")
    done
    echo "All program files:"
    echo "${PROGRAM_FILES[@]}"
    export PROGRAM_FILES
fi

move_output_files() {
    if [ "$TYPE" = "DMD" ]; then
        echo "${EX_DMD_PATH_LOCAL}"
        echo "${SCRIPT_NAME}"
        cd ${EX_DMD_PATH_LOCAL}
        mkdir -p "${SCRIPT_NAME}_out"
        ALL_FILES=(*)
        for file in "${ALL_FILES[@]}"; do
            movable=true
            for p in "${PROGRAM_FILES[@]}"; do
                f=($(basename $file))
                if [[ "$f" = "$p" || "$f" = *_out* ]]; then
                    movable=false
                    break
                fi
            done
            if [[ $movable = true ]]; then
                mv "${f}" "${EX_DMD_PATH_LOCAL}/${SCRIPT_NAME}_out"
            fi
        done

        cd ${EX_DMD_PATH_BASELINE}
        mkdir -p "${SCRIPT_NAME}_out"
        ALL_FILES=(*)
        for file in "${ALL_FILES[@]}"; do
            movable=true
            for p in "${PROGRAM_FILES[@]}"; do
                f=($(basename $file))
                if [[ "$f" = "$p" || "$f" = *_out* ]]; then
                    movable=false
                    break
                fi
            done
            if [[ "$movable" = true ]]; then
                mv "${f}" "${EX_DMD_PATH_BASELINE}/${SCRIPT_NAME}_out"
            fi
        done
    elif [ "$TYPE" = "PROM" ]; then
        cd ${EX_PROM_PATH_LOCAL}
        mkdir -p "${SCRIPT_NAME}_out"
        ALL_FILES=(*)
        for file in "${ALL_FILES[@]}"; do
            movable=true
            for p in "${PROGRAM_FILES[@]}"; do
                f=($(basename $file))
                if [[ "$f" = "$p" || "$f" = *_out* ]]; then
                    movable=false
                    break
                fi
            done
            if [[ $movable = true ]]; then
                mv "${f}" "${EX_PROM_PATH_LOCAL}/${SCRIPT_NAME}_out"
            fi
        done

        echo "${EX_PROM_PATH_BASELINE}"
        echo "${SCRIPT_NAME}"
        cd ${EX_DMD_PATH_BASELINE}
        mkdir -p "${SCRIPT_NAME}_out"
        ALL_FILES=(*)
        for file in "${ALL_FILES[@]}"; do
            movable=true
            for p in "${PROGRAM_FILES[@]}"; do
                f=($(basename $file))
                if [[ "$f" = "$p" || "$f" = *_out* ]]; then
                    movable=false
                    break
                fi
            done
            if [[ "$movable" = true ]]; then
                mv "${f}" "${EX_PROM_PATH_BASELINE}/${SCRIPT_NAME}_out"
            fi
        done
    else echo "ERROR: neither DMD nor PROM"
    fi
}

run_cmds() {
    for cmd in "${CMDS[@]}"; do
        eval "$cmd"
    done
}
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
