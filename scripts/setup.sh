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

export CFLAGS="-fPIC ${CFLAGS}"
export CPPFLAGS="-fPIC ${CPPFLAGS}"
export CXXFLAGS="-fPIC ${CXXFLAGS}"

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
  # NOTE(kevin): This is the hypre version used by PyMFEM v4.5.2.0.
  # This version matching is required to support pylibROM-PyMFEM interface.
  wget https://github.com/hypre-space/hypre/archive/v2.28.0.tar.gz
  tar -zxvf v2.28.0.tar.gz
  mv hypre-2.28.0 hypre
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
  make config shared=1
  check_result $? parmetis-config
  make -j 8
  check_result $? parmetis-installation
  METIS_DIR=$LIB_DIR/parmetis-4.0.3
  METIS_OPT=-I${METIS_DIR}/metis/include
  cd ${METIS_DIR}/build
  MACHINE_ARCH=$(ls)
  ln -s $MACHINE_ARCH lib
  check_result $? parmetis-link
fi


METIS_DIR=$LIB_DIR/parmetis-4.0.3
METIS_OPT=-I${METIS_DIR}/metis/include
METIS_LIB="-L${METIS_DIR}/build/lib/libparmetis -lparmetis -L${METIS_DIR}/build/lib/libmetis -lmetis"

if [ ${MFEM_USE_LAPACK} == "On" ]; then
  # MFEM makefile needs YES/NO
  MFEM_USE_LAPACK="YES"
fi

# Install distributed version of SuperLU
cd $LIB_DIR
if [ ! -z ${INSTALL_SUPERLU} ] && [ ! -d "superlu_dist" ]; then
  wget https://github.com/xiaoyeli/superlu_dist/archive/refs/tags/v6.3.1.tar.gz
  tar -zxvf v6.3.1.tar.gz
  mv superlu_dist-6.3.1 superlu_dist
  cd superlu_dist
  mkdir build
  cd build
  export PARMETIS_ROOT=${METIS_DIR}
  export PARMETIS_BUILD_DIR=${PARMETIS_ROOT}/build/lib/
  cmake -DCMAKE_INSTALL_PREFIX=./ \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_C_COMPILER=mpicc \
        -DCMAKE_CXX_COMPILER=mpicxx \
        -DCMAKE_Fortran_COMPILER=mpifort \
        -Denable_complex16=OFF \
        -Denable_examples=OFF \
        -DTPL_ENABLE_PARMETISLIB=ON \
        -DTPL_PARMETIS_INCLUDE_DIRS="${PARMETIS_ROOT}/include;${PARMETIS_ROOT}/metis/include" \
        -DTPL_PARMETIS_LIBRARIES="${PARMETIS_BUILD_DIR}/libparmetis/libparmetis.so;${PARMETIS_BUILD_DIR}/libmetis/libmetis.so" \
        -DTPL_ENABLE_COMBBLASLIB=OFF \
        -DUSE_XSDK_DEFAULTS="False" \
        -DBUILD_SHARED_LIBS=ON \
        -DBUILD_STATIC_LIBS=OFF \
        ..
    check_result $? superlu-cmake
    make -j 8
    check_result $? superlu-build
    make install
    check_result $? superlu-installation
fi

SUPERLU_DIR=$LIB_DIR/superlu_dist
SUPERLU_OPT=-I${SUPERLU_DIR}/build/include
SUPERLU_LIB="-Wl,-rpath,${SUPERLU_DIR}/build/lib/ -L${SUPERLU_DIR}/build/lib/ -lsuperlu_dist -lblas"

# backup any user provided flags and clear them before compiling MFEM
ROM_CFLAGS="${CFLAGS}"
ROM_CXXFLAGS="${CXXFLAGS}"
unset CFLAGS
unset CXXFLAGS

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
        make -j 8 pdebug CPPFLAGS="${ROM_CXXFLAGS}" STATIC=NO SHARED=YES MFEM_USE_MPI=YES MFEM_USE_GSLIB=${MG} MFEM_USE_LAPACK=${MFEM_USE_LAPACK:-"NO"} MFEM_USE_METIS=YES MFEM_USE_METIS_5=YES METIS_DIR="$METIS_DIR" METIS_OPT="$METIS_OPT" METIS_LIB="$METIS_LIB" MFEM_USE_SUPERLU=${MFEM_USE_SUPERLU:-"NO"} SUPERLU_DIR="$SUPERLU_DIR" SUPERLU_OPT="$SUPERLU_OPT" SUPERLU_LIB="$SUPERLU_LIB"
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
        # v4.7.0 commit. This is the mfem version used by PyMFEM v4.7.0.1.
        # This version matching is required to support pylibROM-PyMFEM interface.
        git checkout dc9128ef596e84daf1138aa3046b826bba9d259f
        make -j 8 parallel CPPFLAGS="${ROM_CXXFLAGS}" STATIC=NO SHARED=YES MFEM_USE_MPI=YES MFEM_USE_GSLIB=${MG} MFEM_USE_LAPACK=${MFEM_USE_LAPACK:-"NO"} MFEM_USE_METIS=YES MFEM_USE_METIS_5=YES METIS_DIR="$METIS_DIR" METIS_OPT="$METIS_OPT" METIS_LIB="$METIS_LIB" MFEM_USE_SUPERLU=${MFEM_USE_SUPERLU:-"NO"} SUPERLU_DIR="$SUPERLU_DIR" SUPERLU_OPT="$SUPERLU_OPT" SUPERLU_LIB="$SUPERLU_LIB"
        check_result $? mfem-parallel-installation
    fi
    cd $LIB_DIR
    rm mfem
    ln -s mfem_parallel mfem
    check_result $? mfem-parallel-link
fi

# restore any user provided flags for compiling libROM
export CFLAGS="${ROM_CFLAGS}"
export CXXFLAGS="${ROM_CXXFLAGS}"
