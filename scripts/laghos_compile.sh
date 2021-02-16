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
pushd ${REPO_PREFIX}/build
rm -rf *
mkdir -p ${REPO_PREFIX}/dependencies

if [ "$(uname)" == "Darwin" ]; then
  which -s brew > /dev/null
  if [[ $? != 0 ]] ; then
      # Install Homebrew
      echo "Homebrew installation is required."
      exit 1
  fi
  xcode-select -p > /dev/null
  if [[ $? != 0 ]] ; then
      xcode-select --install
  fi
  brew list open-mpi > /dev/null || brew install open-mpi
  brew list openblas > /dev/null || brew install openblas
  brew list lapack > /dev/null || brew install lapack
  brew list scalapack > /dev/null || brew install scalapack
  brew list hdf5 > /dev/null || brew install hdf5
  if [[ $1 == "--torch" ]];
  then
    if [[ ! -d ${REPO_PREFIX}/dependencies/libtorch ]];
    then
      wget -P ${REPO_PREFIX}/dependencies https://download.pytorch.org/libtorch/cpu/libtorch-macos-1.6.0.zip
      unzip ${REPO_PREFIX}/dependencies/libtorch-macos-1.6.0.zip -d ${REPO_PREFIX}/dependencies
      install_name_tool -id @rpath/libiomp5.dylib ${REPO_PREFIX}/dependencies/libtorch/lib/libiomp5.dylib
    fi
  fi
  cmake ${REPO_PREFIX}
  make
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
  TOOLCHAIN_FILE=${REPO_PREFIX}/cmake/toolchains/ic19-toss_3_x86_64_ib-librom-dev.cmake
  if [[ $1 == "--torch" ]];
  then
    if [[ ! -d ${REPO_PREFIX}/dependencies/libtorch ]];
    then
      wget -P ${REPO_PREFIX}/dependencies https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-1.6.0%2Bcpu.zip
      unzip ${REPO_PREFIX}/dependencies/libtorch-shared-with-deps-1.6.0+cpu.zip -d ${REPO_PREFIX}/dependencies
    fi
  fi
  cmake ${REPO_PREFIX} \
        -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN_FILE} \
        -DCMAKE_BUILD_TYPE=Optimized
  make VERBOSE=1 -j8
fi
popd
