#!/usr/bin/env bash

export PKG_ROOT=/usr/gapps/bdiv/toss_3_x86_64_ib/clang-3.9-mvapich2-2.2
export COMPILER_ROOT=/usr/tce/packages/mvapich2/mvapich2-2.2-clang-3.9.1/bin
export LAPACK_LDFLAGS="-L${PKG_ROOT}/lapack/3.5.0/lib -llapack -lblas"
./configure \
    --with-CXX=${COMPILER_ROOT}/mpicxx \
    --with-FC=${COMPILER_ROOT}/mpif90 \
    --with-lapack=${PKG_ROOT}/lapack/3.5.0 \
    --with-lapack-libs="${LAPACK_LDFLAGS}" \
    --with-hdf5=${PKG_ROOT}/hdf5/1.8.10p1 \
    --with-zlib=${PKG_ROOT}/zlib/1.2.3/lib \
    --with-gtest=${PKG_ROOT}/gtest/1.8.0 \
    --with-elemental=no \
    --enable-opt=yes \
    --enable-debug=no \
    "$@"
