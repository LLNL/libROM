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

#NOTE: Spectrum MPI is IBM's prepackaged version of OpenMPI
#NOTE: Extra LD flags are what is used by /usr/tcetmp/bin/xlflang -v
#      when called on a Fortran 90 source file. ALCF recommends
#      including these flags when symbols like _xldipow are not found
#      See https://wiki.alcf.anl.gov/old/index.php/Compiling_and_linking

export PKG_ROOT=/usr/gapps/bdiv/blueos_3_ppc64le_ib/clang-6-openmpi
export SPECTRUM_PATH=/usr/tce/packages/spectrum-mpi
export SPECTRUM_VERSION=spectrum-mpi-rolling-release-xl-beta-2018.08.08
export COMPILER_ROOT=${SPECTRUM_PATH}/${SPECTRUM_VERSION}/bin

# Either of these pairs of BLAS/LAPACK libraries will work

#export LAPACK_ROOT=${PKG_ROOT}/lapack/3.5.0
#export BLAS_ROOT=${PKG_ROOT}/lapack/3.5.0
export LAPACK_ROOT=/usr/tcetmp/packages/lapack/lapack-3.6.0-xlf-15.1.5
export BLAS_ROOT=/usr/tcetmp/packages/blas/blas-3.6.0-xlf-15.1.5

# Use full paths to force static linkage; the usual -L/-l flags default to
# linking shared libraries, and the beta version of xlf seems like it
# doesn't really support RPATH. Whether or not xlf supports RPATH,
# the usual -L/-l flags don't work for LAPACK/BLAS
export LAPACK_LDFLAGS="${LAPACK_ROOT}/lib/liblapack.a ${BLAS_ROOT}/lib/libblas.a"

# Use xlflang because Fortran packages won't compile if mpixlf used
export CC=${COMPILER_ROOT}/mpicc
./configure \
    --with-CXX=${COMPILER_ROOT}/mpicxx \
    --with-FC="/usr/tce/packages/xl/xl-beta-2018.08.08/bin/xlf90" \
    --with-lapack=${LAPACK_ROOT} \
    --with-lapack-libs="${LAPACK_LDFLAGS}" \
    --with-hdf5=${PKG_ROOT}/hdf5/1.8.10p1 \
    --with-zlib=${PKG_ROOT}/zlib/1.2.3/lib \
    --with-gtest=no \
    --with-extra-ld-flags="-L/usr/tce/packages/xl/xl-beta-2018.08.08/xlf/16.1.1/lib -lxl -lxlf90 -lxlf90_r -lxlfmath -lm -lxlopt" \
    --with-elemental=no \
    --enable-opt=yes \
    --enable-debug=no \
    "$@"
