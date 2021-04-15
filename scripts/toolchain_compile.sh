#!/usr/bin/env bash

###############################################################################
#
#  Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
#  and other libROM project developers. See the top-level COPYRIGHT
#  file for details.
#
#  SPDX-License-Identifier: (Apache-2.0 OR MIT)
#
###############################################################################

REPO_PREFIX=$(git rev-parse --show-toplevel)

if [ "$1" == "" ] || [ $# -gt 1 ]; then
    echo "Usage: ./toolchain_compile.sh toolchain.cmake"
    echo "toolchain.cmake should be placed in /path/to/librom/cmake/toolchains"
    exit 0
fi

TOOLCHAIN_FILE=${REPO_PREFIX}/cmake/toolchains/$1
pushd ${REPO_PREFIX}/build
rm -rf *

cmake ${REPO_PREFIX} \
      -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN_FILE} \
      -DCMAKE_BUILD_TYPE=Optimized
make VERBOSE=1 -j8
popd
