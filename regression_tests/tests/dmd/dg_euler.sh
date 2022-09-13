#!/bin/bash
source $GITHUB_WORKSPACE/regression_tests/common.sh

cd ${GITHUB_WORKSPACE}/build/examples/dmd
mpirun -np 8 dg_euler -p 1 -rs 1 -rp 1 -o 5 -s 6 -tf 0.1

cd ${BASELINE_DIR}/libROM/build/examples/dmd # Baseline(master) branch libROM
mpirun -np 8 dg_euler -p 1 -rs 1 -rp 1 -o 5 -s 6 -tf 0.1

cd ${GITHUB_WORKSPACE}/build/tests

echo "Running solution comparator for $0 with $1 processors"
#./basisComparator ${GITHUB_WORKSPACE}/build/examples/dmd/dg_euler ${BASELINE_DIR}/libROM/build/examples/dmd/dg_euler "1e-7" "$1"
./solutionComparator "${GITHUB_WORKSPACE}/build/examples/dmd/vortex-0-final.000000" "${BASELINE_DIR}/libROM/build/examples/dmd/vortex-0-final.000000" "1.0e-5" "$1"
./solutionComparator "${GITHUB_WORKSPACE}/build/examples/dmd/vortex-1-final.000000" "${BASELINE_DIR}/libROM/build/examples/dmd/vortex-1-final.000000" "1.0e-5" "$1"
./solutionComparator "${GITHUB_WORKSPACE}/build/examples/dmd/vortex-2-final.000000" "${BASELINE_DIR}/libROM/build/examples/dmd/vortex-2-final.000000" "1.0e-5" "$1"
./solutionComparator "${GITHUB_WORKSPACE}/build/examples/dmd/vortex-3-final.000000" "${BASELINE_DIR}/libROM/build/examples/dmd/vortex-3-final.000000" "1.0e-5" "$1"












