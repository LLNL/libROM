#!/usr/bin/env bash

export PKG_ROOT=/usr/gapps/bdiv/toss_3_x86_64_ib/intel-18-mvapich2-2.2
export COMPILER_ROOT=/usr/tce/packages/mvapich2/mvapich2-2.2-intel-18.0.2/bin
export GCC_ROOT=/usr/tce/packages/gcc/gcc-7.1.0/bin
./configure \
    --with-CXX="${COMPILER_ROOT}/mpicxx -gxx-name=${GCC_ROOT}/g++" \
    --with-FC="${COMPILER_ROOT}/mpif90" \
    --with-lapack=${PKG_ROOT}/lapack/3.5.0 \
    --with-hdf5=${PKG_ROOT}/hdf5/1.8.10p1 \
    --with-zlib=${PKG_ROOT}/zlib/1.2.3/lib \
    --with-gtest=${PKG_ROOT}/gmock/1.7.0/gtest \
    --with-elemental=no \
    --enable-opt=yes \
    --enable-debug=no \
    "$@"
