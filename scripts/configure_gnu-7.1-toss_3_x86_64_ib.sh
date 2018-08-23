#!/usr/bin/env bash

export PKG_ROOT=/usr/gapps/bdiv/toss_3_x86_64_ib/gnu-7.1-mvapich2-2.2
export COMPILER_ROOT=/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-7.1.0/bin
export LAPACK_LDFLAGS="-L${PKG_ROOT}/lapack/3.5.0/lib -llapack -lblas"
# Needed to fix failing CPP check
#export CC="${COMPILER_ROOT}/mpicc -gxx-name=${GCC_ROOT}/g++"

./configure \
    --with-FC="${COMPILER_ROOT}/mpif90" \
    --with-CXX="${COMPILER_ROOT}/mpicxx" \
    --with-lapack=${PKG_ROOT}/lapack/3.5.0/lib \
    --with-lapack-libs="${LAPACK_LDFLAGS}" \
    --with-hdf5=${PKG_ROOT}/hdf5/1.8.10p1 \
    --with-zlib=${PKG_ROOT}/zlib/1.2.3/lib \
    --with-gtest=no \
    --with-elemental=no \
    --enable-opt=yes \
    --enable-debug=no \
    "$@"
