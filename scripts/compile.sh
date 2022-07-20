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
UPDATE_LIBS=false
INSTALL_PYTHON=false

# Get options
while getopts "ah:dh:mh:t:uh:ph" o;
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
            INSTALL_PYTHON=true
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
    . ${REPO_PREFIX}/scripts/mfem_setup.sh
fi

if [[ $ARDRA == "true" ]]; then
    mkdir -p ${REPO_PREFIX}/buildArdra
    LIB_BUILD_DIR=${REPO_PREFIX}/buildArdra/lib
    pushd ${REPO_PREFIX}/buildArdra
else
    mkdir -p ${REPO_PREFIX}/build
    LIB_BUILD_DIR=${REPO_PREFIX}/build/lib
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
        -DUSE_MFEM=${USE_MFEM}
  make
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
  if [[ $ARDRA == "true" ]]; then
      TOOLCHAIN_FILE=${REPO_PREFIX}/cmake/toolchains/ic18-toss_3_x86_64_ib-ardra.cmake
  elif [[ -z ${TOOLCHAIN_FILE} ]]; then
      TOOLCHAIN_FILE=${REPO_PREFIX}/cmake/toolchains/default-toss_3_x86_64_ib-librom-dev.cmake
  fi
  cmake ${REPO_PREFIX} \
        -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN_FILE} \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -DUSE_MFEM=${USE_MFEM}
  make -j8
fi

popd

if [[ $INSTALL_PYTHON == "true" ]]; then
    . ${REPO_PREFIX}/scripts/python_setup.sh
    cd ${REPO_PREFIX}/lib
    ${REPO_PREFIX}/dependencies/swig/swig_install/bin/swig -I${REPO_PREFIX}/dependencies/swig/swig_install/share/swig/4.0.2 -I${REPO_PREFIX}/dependencies/swig/swig_install/share/swig/4.0.2/python -I${REPO_PREFIX}/dependencies/swig/swig_install/share/swig/4.0.2/std -c++ -python -py3 librom.i
    mv librom_wrap.cxx $LIB_BUILD_DIR
    mv librom.py $LIB_BUILD_DIR
    cd $LIB_BUILD_DIR

    # TODO: Untie the Python lib from LC.
    mpicxx -O2 -fPIC -DSWIG -c librom_wrap.cxx -I/usr/include/python3.6m -I/usr/workspace/huynh24/libROM4/env3/lib/python3.7/site-packages/mpi4py/include -I${REPO_PREFIX}/lib
    mpicxx -shared -L/usr/workspace/huynh24/libROM4/env3/lib/python3.7/site-packages -Wl,-rpath,$LIB_BUILD_DIR -L$LIB_BUILD_DIR librom_wrap.o -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl -lROM -o _librom.so
fi
