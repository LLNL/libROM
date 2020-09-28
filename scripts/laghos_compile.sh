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
  cmake ${REPO_PREFIX}
  make
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
  TOOLCHAIN_FILE=${REPO_PREFIX}/cmake/toolchains/ic19-toss_3_x86_64_ib-librom-dev.cmake
  cmake ${REPO_PREFIX} \
        -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN_FILE} \
        -DCMAKE_BUILD_TYPE=Optimized \
        "$@"
  make VERBOSE=1 -j8
fi
popd
