#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

CMDS=( 
    "$COMMAND dg_euler -p 1 -rs 1 -rp 1 -o 5 -s 6 -tf 0.1" 
)
TYPE="DMD"
cd ${EX_DMD_PATH_LOCAL}
run_cmds

cd ${EX_DMD_PATH_BASELINE} # Baseline(master) branch libROM
run_cmds

cd ${GITHUB_WORKSPACE}/build/tests

echo "Running solution comparator for $0 with $1 processors"
#./basisComparator ${GITHUB_WORKSPACE}/build/examples/dmd/dg_euler ${BASELINE_DIR}/libROM/build/examples/dmd/dg_euler "1e-7" "$1"

./solutionComparator "${EX_DMD_PATH_LOCAL}/vortex-0-final.000000" "${EX_DMD_PATH_BASELINE}/vortex-0-final.000000" "1.0e-5" "8"
check_fail

./solutionComparator "${EX_DMD_PATH_LOCAL}/vortex-1-final.000000" "${EX_DMD_PATH_BASELINE}/vortex-1-final.000000" "1.0e-5" "$1"
check_fail

./solutionComparator "${EX_DMD_PATH_LOCAL}/vortex-2-final.000000" "${EX_DMD_PATH_BASELINE}/vortex-2-final.000000" "1.0e-5" "$1"
check_fail

./solutionComparator "${EX_DMD_PATH_LOCAL}/vortex-3-final.000000" "${EX_DMD_PATH_BASELINE}/vortex-3-final.000000" "1.0e-5" "$1"
check_fail

move_output_files












