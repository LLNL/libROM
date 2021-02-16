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
TOOLCHAIN_FILE=${REPO_PREFIX}/cmake/toolchains/ic18-toss_3_x86_64_ib-ardra.cmake

mkdir ${REPO_PREFIX}/buildArdra
pushd ${REPO_PREFIX}/buildArdra
rm -rf *
mkdir -p ${REPO_PREFIX}/dependencies
if [[ $1 == "--torch" ]];
then
  if [[ ! -d ${REPO_PREFIX}/libtorch ]];
  then
    wget -P ${REPO_PREFIX} https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-1.6.0%2Bcpu.zip
    unzip ${REPO_PREFIX}/libtorch-shared-with-deps-1.6.0+cpu.zip -d ${REPO_PREFIX}
  fi
fi
cmake ${REPO_PREFIX} \
      -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN_FILE} \
      -DCMAKE_BUILD_TYPE=Optimized \
      "$@"
make VERBOSE=1 -j8
popd
