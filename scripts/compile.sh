#!/usr/bin/env bash

###############################################################################
#
#  Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
#  and other libROM project developers. See the top-level COPYRIGHT
#  file for details.
#
#  SPDX-License-Identifier: (Apache-2.0 OR MIT)
#
###############################################################################
check_result () {
  # $1: Result output of the previous command ($?)
  # $2: Name of the previous command
  if [ $1 -eq 0 ]; then
      echo "$2 succeeded"
  else
      echo "$2 failed"
      exit -1
  fi
}

ARDRA=false
BUILD_TYPE="Optimized"
USE_MFEM="Off"
UPDATE_LIBS=false
INSTALL_SCALAPACK=false
MFEM_USE_GSLIB="Off"

cleanup_dependencies() {
    pushd .
    cd "$HOME_DIR"
    rm ./build/CMakeCache.txt
    pushd dependencies || exit 1
    mydirs=(*)
    for dir in "${mydirs[@]}" ; do
        if [[ $dir =~ "mfem" ]]; then
            echo "Deleting $dir"
            rm -rf $dir
        fi
    done
    popd
    popd
}


# Get options
while getopts "ah:dh:gh:mh:t:uh:sh" o;
do
    case "${o}" in
        a)
            ARDRA=true
            ;;
        d)
            BUILD_TYPE="Debug"
            ;;
        g)
            MFEM_USE_GSLIB="On"
            ;;
        m)
            USE_MFEM="On"
            ;;
        t)
            TOOLCHAIN_FILE=${OPTARG}
            ;;
        u)
            UPDATE_LIBS=true
            ;;
        s)
            INSTALL_SCALAPACK=true
            ;;
    *)
            echo "Unknown option."
            exit 1
      ;;
    esac
done
shift $((OPTIND-1))

CURR_DIR=$(pwd) 
if [[ -d "dependencies" ]]; then
    HOME_DIR=$CURR_DIR
elif [[ -f "compile.sh" ]]; then
    HOME_DIR=..
else
    echo "You can only run compile.sh from either the libROM or the scripts directory"
    exit
fi
# If both include and exclude are set, fail
if [[ -n "${TOOLCHAIN_FILE}" ]] && [[ $ARDRA == "true" ]]; then
     echo "Choose only Ardra or add your own toolchain file, not both."
	 exit 1
fi

if [[ $MFEM_USE_GSLIB == "On" ]] && [[ $USE_MFEM == "Off" ]]; then
    echo "gslib can only be used with mfem"
	exit 1
fi

if [[ $MFEM_USE_GSLIB == "On" ]] && [[ ! -d "$HOME_DIR/dependencies/gslib" ]]; then
   cleanup_dependencies
fi
export MFEM_USE_GSLIB

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
REPO_PREFIX=$( dirname $SCRIPT_DIR )

if [[ $USE_MFEM == "On" ]]; then
    . ${REPO_PREFIX}/scripts/setup.sh ${INSTALL_SCALAPACK}
    check_result $? scripts-setup
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
  brew list open-mpi > /dev/null || brew install open-mpi
  brew list openblas > /dev/null || brew install openblas
  brew list lapack > /dev/null || brew install lapack
  brew list scalapack > /dev/null || brew install scalapack
  brew list hdf5 > /dev/null || brew install hdf5
  brew list cmake > /dev/null || brew install cmake
  cmake ${REPO_PREFIX} \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -DUSE_MFEM=${USE_MFEM} \
        -DMFEM_USE_GSLIB=${MFEM_USE_GSLIB}
  check_result $? librom-config
  make
  check_result $? librom-build
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
  if [[ $ARDRA == "true" ]]; then
      TOOLCHAIN_FILE=${REPO_PREFIX}/cmake/toolchains/ic18-toss_3_x86_64_ib-ardra.cmake
  elif [[ -z ${TOOLCHAIN_FILE} ]]; then
      TOOLCHAIN_FILE=${REPO_PREFIX}/cmake/toolchains/default-toss_4_x86_64_ib-librom-dev.cmake
  fi
  cmake ${REPO_PREFIX} \
        -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN_FILE} \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -DUSE_MFEM=${USE_MFEM} \
        -DMFEM_USE_GSLIB=${MFEM_USE_GSLIB}
  check_result $? librom-config
  make -j8
  check_result $? librom-build
fi
popd
