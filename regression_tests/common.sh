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
export OFFSET # Number of header lines
trap "move_output_files_after_error" EXIT
SCRIPT_NAME=$(basename "$0" ".sh")
NP=$(($1))
NUM_PROCESSES=4
# Check machine
#echo "OS: $MACHINE"
if [[ $MACHINE = "Linux" ]]; then
	COMMAND="srun -p pbatch -n $NUM_PROCESSES"
elif [[ $MACHINE = "Darwin" ]]; then
	COMMAND="mpirun -np $NUM_PROCESSES"
elif [[ $MACHINE = "GitHub" ]]; then
    NUM_PROCESSES=1
	COMMAND=""
else
    echo "Bad OS: $MACHINE"
    exit 1
fi

move_output_files_after_error() {
    if [[ "$?" -eq 0 ]]; then
        return
    fi
    move_output_files
    exit 1
}

move_output_files() { 
    EX_PATHS=("${EX_DMD_PATH_LOCAL}" "${EX_DMD_PATH_BASELINE}" "${EX_PROM_PATH_LOCAL}" "${EX_PROM_PATH_BASELINE}")
    for path in "${EX_PATHS[@]}"; do
        if [[ "$TYPE" == "DMD" && ($path ==  "${EX_PROM_PATH_LOCAL}" || $path == "${EX_PROM_PATH_BASELINE}") ]]; then
            continue
        elif [[ "$TYPE" = "PROM" && ($path ==  "${EX_DMD_PATH_LOCAL}" || $path == "${EX_DMD_PATH_BASELINE}") ]]; then
            continue
        fi
        cd "${path}"
        rm -rf "${SCRIPT_NAME}_out"
        mkdir "${SCRIPT_NAME}_out"
        $(find . -maxdepth 1  -not -name "*_out" -not -name "." ! -perm -700 -exec mv {} "${SCRIPT_NAME}_out" \;)
        $(find . -maxdepth 1  -type d -not -name "*_out" -not -name "." -exec mv {} "${SCRIPT_NAME}_out" \;)
    done
}

run_cmds() {
    for cmd in "${CMDS[@]}"; do
        eval "$cmd"
    done
}

run_tests() {
    # Run commands from the local directory
    if [[ $TYPE == "DMD" ]]; then
        cd  ${EX_DMD_PATH_LOCAL}
    elif [[ $TYPE == "PROM" ]]; then
        cd ${EX_PROM_PATH_LOCAL}
    else
        echo "Unrecognized TYPE is ${TYPE}"
    fi
    run_cmds

    # Run commands from the baseline directory
    if [[ $TYPE == "DMD" ]]; then
        cd  ${EX_DMD_PATH_BASELINE}
    elif [[ $TYPE == "PROM" ]]; then
        cd ${EX_PROM_PATH_BASELINE}
    else
        echo "Unrecognized TYPE is ${TYPE}"
    fi
    run_cmds

    # Compare results between the two
    files_to_compare=(*)
    # if the test has FILES defined, compare those files instead
    if [[ ! -z "${FILES}" ]]; then
        files_to_compare=(${FILES[@]})
    fi
    cd ${GITHUB_WORKSPACE}/build/tests

    for f in "${files_to_compare[@]}"; do
        if [[ $f =~ final || "$f" == "sol"*".000000" && "$f" != "sol_dofs"* || "$f" == "Sol0"  ]]; then
            if [[ $TYPE == "DMD" && "$f" == *".0000"* && $MACHINE = "GitHub" ]]; then
                if [[ $OFFSET -ne 0 ]]; then
                    cp "${EX_DMD_PATH_BASELINE}/${f}" "${EX_DMD_PATH_BASELINE}/${f}-orig"
                    cp "${EX_DMD_PATH_LOCAL}/${f}" "${EX_DMD_PATH_LOCAL}/${f}-orig"
                    sed -i  '1,'"$OFFSET"'d' "${EX_DMD_PATH_BASELINE}/${f}"
                    sed -i '1,'"$OFFSET"'d' "${EX_DMD_PATH_LOCAL}/${f}"
                fi          
            elif [[ $TYPE == "DMD" && "$f" == *".0000"* ]]; then
                if [[ $OFFSET -ne 0 ]]; then
                    cp "${EX_DMD_PATH_BASELINE}/${f}" "${EX_DMD_PATH_BASELINE}/${f}-orig"
                    cp "${EX_DMD_PATH_LOCAL}/${f}" "${EX_DMD_PATH_LOCAL}/${f}-orig"
                    if [[ $MACHINE = "Darwin" ]]; then
                        sed -i '' '1,'"$OFFSET"'d' "${EX_DMD_PATH_BASELINE}/${f}"
                        sed -i '' '1,'"$OFFSET"'d' "${EX_DMD_PATH_LOCAL}/${f}"
                    elif [[ $MACHINE == "Linux"  || $MACHINE == "GitHub" ]]; then
                        sed -i '1,'"$OFFSET"'d' "${EX_DMD_PATH_BASELINE}/${f}"
                        sed -i  '1,'"$OFFSET"'d' "${EX_DMD_PATH_LOCAL}/${f}"
                    fi
                fi
            elif [[ $TYPE == "PROM" ]]; then
                if [[ $OFFSET -ne 0 ]]; then
                    cp "${EX_PROM_PATH_BASELINE}/${f}" "${EX_PROM_PATH_BASELINE}/${f}-orig"
                    cp "${EX_PROM_PATH_LOCAL}/${f}" "${EX_PROM_PATH_LOCAL}/${f}-orig"
                    if [[ $MACHINE = "Darwin" ]]; then
                        sed -i '' '1,'"$OFFSET"'d' "${EX_PROM_PATH_BASELINE}/${f}"
                        sed -i '' '1,'"$OFFSET"'d' "${EX_PROM_PATH_LOCAL}/${f}"
                    elif [[ $MACHINE == "Linux" || $MACHINE == "GitHub" ]]; then
                        sed -i '1,'"$OFFSET"'d' "${EX_PROM_PATH_BASELINE}/${f}"
                        sed -i '1,'"$OFFSET"'d' "${EX_PROM_PATH_LOCAL}/${f}"
                    fi
                fi
            else
                continue
            fi
        fi
    done

    for f in "${files_to_compare[@]}"; do
        if [[ ! $f =~ snapshot && $f =~ basis && -n $test_offline ]]; then # Do not compare offline results(bases) by default. Do not compare sampled snapshots
            fn="${f%.*}"
            if [[ $TYPE == "DMD" ]]; then
                ./basisComparator "${EX_DMD_PATH_BASELINE}/$fn" "${EX_DMD_PATH_LOCAL}/$fn" 1e-7 1
            elif [[ $TYPE == "PROM" ]]; then
                ./basisComparator "${EX_PROM_PATH_BASELINE}/$fn" "${EX_PROM_PATH_LOCAL}/$fn" 1e-7 1
            else
                continue
            fi
            check_fail 
        elif [[ $f =~ final || "$f" == "sol"*".000000" && "$f" != "sol_dofs"* || "$f" == "Sol0"  ]]; then
            if [[ $TYPE == "DMD" && "$f" == *".000000" && $MACHINE = "GitHub" ]]; then
                ./solutionComparator "${EX_DMD_PATH_BASELINE}/${f}"  "${EX_DMD_PATH_LOCAL}/${f}" "1.0e-5" "$NUM_PROCESSES"            
            elif [[ $TYPE == "DMD" && "$f" == *".000000" ]]; then
                ./solutionComparator  "${EX_DMD_PATH_BASELINE}/${f}" "${EX_DMD_PATH_LOCAL}/${f}" "1.0e-5" "$NUM_PROCESSES" 
            elif [[ $TYPE == "PROM" ]]; then
                ./solutionComparator "${EX_PROM_PATH_BASELINE}/$f"  "${EX_PROM_PATH_LOCAL}/$f" "1.0e-5" "1" 
            else
                continue
            fi
            check_fail
        elif [[ ! -z "${FILES}" ]]; then
            echo "Comparing custom file ${f}"
            if [[ $TYPE == "DMD" ]]; then
                ./fileComparator "${EX_DMD_PATH_BASELINE}/${f}" "${EX_DMD_PATH_LOCAL}/${f}" "1.0e-5"
            elif [[ $TYPE == "PROM" ]]; then
                ./fileComparator "${EX_PROM_PATH_BASELINE}/${f}" "${EX_PROM_PATH_LOCAL}/${f}" "1.0e-5"
            else
                continue
            fi
            check_fail
        fi
    done
    move_output_files
}

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
