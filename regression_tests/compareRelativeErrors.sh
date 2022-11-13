#!/bin/bash
set -eo pipefail
compareErrors() {
    echo "Comparing relative errors"
    pushd ${GITHUB_WORKSPACE}/regression_tests/results
    line=$(egrep "Relative error | Relative l2 error " ${scriptName}.log)
    if [[ -z $line ]]; then
        echo "Couldn't find any lines containing Relative error in $scriptName"
        exit 1
    fi
    lin=$(egrep "Relative error | Relative l2 error " ${scriptName}.log | awk '{ print $NF }')
    if [[ -z $lin ]]; then
        echo "Couldn't find any lines containing Relative error"
        exit 1
    fi
    num_fields=$(echo $lin | awk '{ print NF }')
    echo "lin = $lin"
    echo "num_fields = $num_fields"
    remainder=$(($num_fields%2))
    echo "remainder = $remainder"
    if [[ $remainder != 0 ]]; then
        echo "Even instances of the Relative error are expected, found $num_fields instances"
        exit 1
    fi
    set_size=$(echo "$num_fields / 2" | bc)
    echo "set_size = $set_size"
    for (( j=1; j <= $set_size; j++ )); do
        err_local=$(echo $lin | awk -v N=$j '{ print $N }')
        err_local=$(echo $err_local|sed 's/e/*10^/g')
        idx=$j+set_size
        err_baseline=$(echo $lin | awk -v N=$idx '{ print $N }')
        err_baseline=$(echo $err_baseline|sed 's/e/*10^/g')
        echo "err_local=$err_local"
        echo "err_baseline=$err_baseline"
        rel_error_tol=1.1
        err_baseline=$(echo "$rel_error_tol*$err_baseline" | bc)
        too_large_error=$(echo "$err_local > $err_baseline" | bc )
        echo "too_large_error = $too_large_error"
        if [[ $too_large_error != 0 ]]; then
            echo "err_local > err_baseline : err_local = $err_local, err_baseline = $err_baseline"
            exit 1
        fi
    done
    popd
}
