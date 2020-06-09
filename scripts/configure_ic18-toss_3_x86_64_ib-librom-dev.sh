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

# export REPO_PREFIX=$(git rev-parse --show-toplevel)
REPO_PREFIX=$(git rev-parse --show-toplevel)
TOOLCHAIN_FILE=${REPO_PREFIX}/cmake/toolchains/ic18-toss_3_x86_64_ib-librom-dev.cmake
# export PKG_ROOT=/usr/tce/packages
# export COMPILER_ROOT=/usr/tce/packages/mvapich2/mvapich2-2.2-intel-18.0.1/bin
# export GCC_ROOT=/usr/tce/packages/gcc/gcc-4.9.3/bin
# export LAPACK_LDFLAGS="-L/usr/lib64 -llapack -lblas"
# Needed to fix failing CPP check
# export CC="${COMPILER_ROOT}/mpicc -gxx-name=${GCC_ROOT}/g++"

# ./configure \
#     --with-FC="${COMPILER_ROOT}/mpif90 -gxx-name=${GCC_ROOT}/g++" \
#     --with-CXX="${COMPILER_ROOT}/mpicxx -gxx-name=${GCC_ROOT}/g++" \
#     --with-lapack=/usr/lib64 \
#     --with-lapack-libs="${LAPACK_LDFLAGS}" \
#     --with-hdf5=${PKG_ROOT}/hdf5/hdf5-serial-1.8.18-intel-18.0.1 \
#     --enable-opt=yes \
#     --enable-debug=no \
#     --with-elemental=no \
#     --with-gtest=/collab/usr/gapps/librom/toss3 \
#     "$@"

pushd ${REPO_PREFIX}/build
rm -rf *
cmake ${REPO_PREFIX} \
      -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN_FILE} \
      "$@"
make VERBOSE=1 -j8
popd
