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

ARDRA=false
BUILD_TYPE="Optimized"
USE_MFEM="Off"
MFEM_DIR=""
HYPRE_DIR=""
UPDATE_LIBS=false

# Get options
while getopts "ah:dh:mh:p:i:t:uh" o;
do
    case "${o}" in
        a)
            ARDRA=true
            ;;
        d)
            BUILD_TYPE="Debug"
            ;;
        m)
            USE_MFEM="On"
            ;;
        p)
            MFEM_DIR=${OPTARG}
            ;;
        i)
            HYPRE_DIR=${OPTARG}
            ;;
        t)
            TOOLCHAIN_FILE=${OPTARG}
            ;;
        u)
            UPDATE_LIBS=true
            ;;
    *)
            echo "Unknown option."
            exit 1
      ;;
    esac
done
shift $((OPTIND-1))

# If both include and exclude are set, fail
if [[ -n "${TOOLCHAIN_FILE}" ]] && [[ $ARDRA == "true" ]]; then
    echo "Choose only Ardra or add your own toolchain file, not both."
		exit 1
fi

REPO_PREFIX=$(git rev-parse --show-toplevel)

if [[ $USE_MFEM == "On" ]]; then
    . ${REPO_PREFIX}/scripts/setup.sh
fi

if [[ $ARDRA == "true" ]]; then
    mkdir -p ${REPO_PREFIX}/buildArdra
    pushd ${REPO_PREFIX}/buildArdra
else
    mkdir -p ${REPO_PREFIX}/build
    pushd ${REPO_PREFIX}/build
fi
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
  softwareupdate --all --install --force
  brew list open-mpi > /dev/null || brew install open-mpi
  brew list openblas > /dev/null || brew install openblas
  brew list lapack > /dev/null || brew install lapack
  brew list scalapack > /dev/null || brew install scalapack
  brew list hdf5 > /dev/null || brew install hdf5
  brew list cmake > /dev/null || brew install cmake
  cmake ${REPO_PREFIX} \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -DUSE_MFEM=${USE_MFEM} \
        -DMFEM_DIR=${MFEM_DIR} \
        -DHYPRE_DIR=${HYPRE_DIR}
  make
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
  if [[ $ARDRA == "true" ]]; then
      TOOLCHAIN_FILE=${REPO_PREFIX}/cmake/toolchains/ic18-toss_3_x86_64_ib-ardra.cmake
  elif [[ -z ${TOOLCHAIN_FILE} ]]; then
      TOOLCHAIN_FILE=${REPO_PREFIX}/cmake/toolchains/ic19-toss_3_x86_64_ib-librom-dev.cmake
  fi
  cmake ${REPO_PREFIX} \
        -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN_FILE} \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -DUSE_MFEM=${USE_MFEM} \
        -DMFEM_DIR=${MFEM_DIR} \
        -DHYPRE_DIR=${HYPRE_DIR}
  make -j8
fi
popd
