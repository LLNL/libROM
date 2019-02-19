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

export PKG_ROOT=/usr/gapps/bdiv/toss_3_x86_64_ib/intel-18-mvapich2-2.2
export COMPILER_ROOT=/usr/tce/packages/mvapich2/mvapich2-2.2-intel-18.0.2/bin
export GCC_ROOT=/usr/tce/packages/gcc/gcc-7.1.0/bin
export LAPACK_LDFLAGS="-Wl,-rpath,${PKG_ROOT}/lapack/3.5.0/lib -L${PKG_ROOT}/lapack/3.5.0/lib -llapack -lblas -lm"
# Needed to fix failing CPP check
export CC="${COMPILER_ROOT}/mpicc -gxx-name=${GCC_ROOT}/g++"

./configure \
    --with-FC="${COMPILER_ROOT}/mpif90 -gxx-name=${GCC_ROOT}/g++" \
    --with-CXX="${COMPILER_ROOT}/mpicxx -gxx-name=${GCC_ROOT}/g++" \
    --with-lapack=${PKG_ROOT}/lapack/3.5.0 \
    --with-lapack-libs="${LAPACK_LDFLAGS}" \
    --with-hdf5=${PKG_ROOT}/hdf5/1.8.10p1 \
    --with-zlib=${PKG_ROOT}/zlib/1.2.3/lib \
    --with-gtest=no \
    --with-elemental=no \
    --enable-opt=yes \
    --enable-debug=no \
    "$@"
