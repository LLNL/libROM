#!/bin/bash
compareErrors() {
    echo "Comparing relative errors" >>$simulationLogFile
    pushd ${GITHUB_WORKSPACE}/regression_tests/results > /dev/null
    if [[ $type_of_test == 'dmd' ]]; then
       lin=$(grep "Relative error " ${scriptName}.log | awk '{ print $NF }')
    elif [[ $type_of_test == 'prom' ]]; then
        lin=$(grep -e "Relative l2 error" -e "Relative error of ROM" -e "The relative error is" ${scriptName}.log | awk '{ print $NF }')
    else
        echo "Neither prom nor dmd" >>$simulationLogFile
        popd > /dev/null
        return 1
    fi
    if [[ -z $lin ]]; then
        echo "Couldn't find any lines containing Relative error in $scriptName" >>$simulationLogFile
        popd > /dev/null
        return 1
    fi
    num_fields=$(echo $lin | awk '{ print NF }')
    remainder=$(($num_fields%2))
    if [[ $remainder != 0 ]]; then
        echo "Even instances of the Relative error are expected, found $num_fields instances" >>$simulationLogFile
        popd > /dev/null
        return 1
    fi
    # Pairwise comparison of the relative errors between the local and baseline runs
    set_size=$(echo "$num_fields / 2" | bc)
    for (( j=1; j <= $set_size; j++ )); do
        err_local=$(echo $lin | awk -v N=$j '{ print $N }')
        # Convert scientific notation to decimal
        err_local=$(echo $err_local|sed 's/e/*10^/g')
        idx=$(( $j+$set_size ))
        err_baseline=$(echo $lin | awk -v N=$idx '{ print $N }')
        err_baseline=$(echo $err_baseline|sed 's/e/*10^/g')
        rel_error_tol=1.1
        err_baseline=$(echo "$rel_error_tol * $err_baseline" | bc -l)
        too_large_error=$(echo "$err_local > $err_baseline" | bc -l )
        if [[ $too_large_error != 0 ]]; then
            echo "Relative error comparison: FAIL" >>$simulationLogFile
            echo "err_local > err_baseline : err_local = $err_local, err_baseline (with 10% tolerance) = $err_baseline" >>$simulationLogFile
            popd > /dev/null
            return 1
        fi
    done
    echo "Relative error comparison: PASS" >>$simulationLogFile
    popd > /dev/null
}
