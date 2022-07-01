#!/bin/bash

# Check whether Homebrew or wget is installed
if [ "$(uname)" == "Darwin" ]; then
  which -s brew > /dev/null
  if [[ $? != 0 ]] ; then
      # Install Homebrew
      echo "Homebrew installation is required."
      exit 1
  fi

  which -s wget > /dev/null
  # Install wget
  if [[ $? != 0 ]] ; then
      brew install wget
  fi
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
LIB_DIR=$SCRIPT_DIR/../dependencies
mkdir -p $LIB_DIR

export CFLAGS="-fPIC"
export CPPFLAGS="-fPIC"
export CXXFLAGS="-fPIC"

# Install HYPRE
cd $LIB_DIR
if [ ! -d "hypre" ]; then
  wget https://github.com/hypre-space/hypre/archive/refs/tags/v2.20.0.tar.gz
  tar -zxvf v2.20.0.tar.gz
  mv hypre-2.20.0 hypre
  cd hypre/src
  ./configure --disable-fortran
  make -j
fi

# Install PARMETIS 4.0.3
cd $LIB_DIR
if [ ! -d "parmetis-4.0.3" ]; then

  wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
  tar -zxvf parmetis-4.0.3.tar.gz
  cd parmetis-4.0.3
  make config
  make
  METIS_DIR=$LIB_DIR/parmetis-4.0.3
  METIS_OPT=-I${METIS_DIR}/metis/include
  cd ${METIS_DIR}/build
  MACHINE_ARCH=$(ls)
  ln -s $MACHINE_ARCH lib
fi

unset CFLAGS
unset CPPFLAGS
unset CXXFLAGS

METIS_DIR=$LIB_DIR/parmetis-4.0.3
METIS_OPT=-I${METIS_DIR}/metis/include
METIS_LIB="-L${METIS_DIR}/build/lib/libparmetis -lparmetis -L${METIS_DIR}/build/lib/libmetis -lmetis"

# Install MFEM
cd $LIB_DIR
if [[ $BUILD_TYPE == "Debug" ]]; then
    if [ ! -d "mfem_debug" ]; then
    UPDATE_LIBS=true
    git clone https://github.com/mfem/mfem.git mfem_debug
    fi
    if [[ $UPDATE_LIBS == "true" ]]; then
        cd mfem_debug
        git pull
        make pdebug -j8 STATIC=NO SHARED=YES MFEM_USE_MPI=YES MFEM_USE_METIS=YES MFEM_USE_METIS_5=YES METIS_DIR="$METIS_DIR" METIS_OPT="$METIS_OPT" METIS_LIB="$METIS_LIB"
    fi
    cd $LIB_DIR
    rm mfem
    ln -s mfem_debug mfem
else
    if [ ! -d "mfem_parallel" ]; then
      UPDATE_LIBS=true
      git clone https://github.com/mfem/mfem.git mfem_parallel
    fi
    if [[ $UPDATE_LIBS == "true" ]]; then
        cd mfem_parallel
        git pull
        make parallel -j8 STATIC=NO SHARED=YES MFEM_USE_MPI=YES MFEM_USE_METIS=YES MFEM_USE_METIS_5=YES METIS_DIR="$METIS_DIR" METIS_OPT="$METIS_OPT" METIS_LIB="$METIS_LIB"
    fi
    cd $LIB_DIR
    rm mfem
    ln -s mfem_parallel mfem
fi
