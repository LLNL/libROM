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
TOOLCHAIN_FILE=${REPO_PREFIX}/cmake/toolchains/ic19-toss_3_x86_64_ib-librom-dev.cmake

pushd ${REPO_PREFIX}/build
rm -rf *
cmake ${REPO_PREFIX} \
      -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN_FILE} \
      -DCMAKE_BUILD_TYPE=Optimized \
      "$@"
make VERBOSE=1 -j8
popd
