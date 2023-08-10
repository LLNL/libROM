#!/bin/bash
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

# Take the first argument as option for installing ScaLAPACK.
INSTALL_SCALAPACK=$1

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

# Install ScaLAPACK if specified.
cd $LIB_DIR
if [[ $INSTALL_SCALAPACK == "true" ]]; then
  if [ -f "scalapack-2.2.0/libscalapack.a" ]; then
    echo "Using dependencies/scalapack-2.2.0/libscalapack.a"
  else
    echo "ScaLAPACK is needed!"
    tar -zxvf scalapack-2.2.0.tar.gz
    cp SLmake.inc scalapack-2.2.0/
    cd scalapack-2.2.0/
    make
    check_result $? ScaLAPACK-installation
  fi
fi

# Install HYPRE
cd $LIB_DIR
if [ ! -d "hypre" ]; then
  wget https://github.com/hypre-space/hypre/archive/refs/tags/v2.20.0.tar.gz
  tar -zxvf v2.20.0.tar.gz
  mv hypre-2.20.0 hypre
  cd hypre/src
  ./configure --disable-fortran
  make -j
  check_result $? hypre-installation
fi

# Install GSLIB
cd $LIB_DIR
if [ $MFEM_USE_GSLIB == "On" ]; then
   MG=YES
else
   MG=NO
fi
if [ $MFEM_USE_GSLIB == "On" ] && [ ! -d "gslib" ]; then
  wget https://github.com/gslib/gslib/archive/v1.0.7.tar.gz
  tar -zxvf v1.0.7.tar.gz
  mv gslib-1.0.7 gslib
  cd gslib
  make CC=mpicc -j
  check_result $? gslib-installation
fi

# Install PARMETIS 4.0.3
cd $LIB_DIR
if [ ! -d "parmetis-4.0.3" ]; then

  tar -zxvf parmetis-4.0.3.tar.gz
  cd parmetis-4.0.3
  make config
  check_result $? parmetis-config
  make
  check_result $? parmetis-installation
  METIS_DIR=$LIB_DIR/parmetis-4.0.3
  METIS_OPT=-I${METIS_DIR}/metis/include
  cd ${METIS_DIR}/build
  MACHINE_ARCH=$(ls)
  ln -s $MACHINE_ARCH lib
  check_result $? parmetis-link
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
        make pdebug -j 8 STATIC=NO SHARED=YES MFEM_USE_MPI=YES MFEM_USE_GSLIB=${MG} MFEM_USE_METIS=YES MFEM_USE_METIS_5=YES METIS_DIR="$METIS_DIR" METIS_OPT="$METIS_OPT" METIS_LIB="$METIS_LIB"
        check_result $? mfem-debug-installation
    fi
    cd $LIB_DIR
    rm mfem
    ln -s mfem_debug mfem
    check_result $? mfem-debug-link
else
    if [ ! -d "mfem_parallel" ]; then
      UPDATE_LIBS=true
      git clone https://github.com/mfem/mfem.git mfem_parallel
    fi
    if [[ $UPDATE_LIBS == "true" ]]; then
        cd mfem_parallel
        git pull
        make parallel -j 8 STATIC=NO SHARED=YES MFEM_USE_MPI=YES MFEM_USE_GSLIB=${MG} MFEM_USE_METIS=YES MFEM_USE_METIS_5=YES METIS_DIR="$METIS_DIR" METIS_OPT="$METIS_OPT" METIS_LIB="$METIS_LIB"
        check_result $? mfem-parallel-installation
    fi
    cd $LIB_DIR
    rm mfem
    ln -s mfem_parallel mfem
    check_result $? mfem-parallel-link
fi
