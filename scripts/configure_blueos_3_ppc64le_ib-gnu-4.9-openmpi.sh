#!/usr/bin/env bash

#NOTE: Spectrum MPI is IBM's prepackaged version of OpenMPI
#NOTE: Extra LD flags are what is used by /usr/tcetmp/bin/xlflang -v
#      when called on a Fortran 90 source file. ALCF recommends
#      including these flags when symbols like _xldipow are not found
#      See https://wiki.alcf.anl.gov/old/index.php/Compiling_and_linking

export PKG_ROOT=/usr/gapps/bdiv/blueos_3_ppc64le_ib/gcc-4.9-openmpi
export SPECTRUM_PATH=/usr/tce/packages/spectrum-mpi
export SPECTRUM_GCC=${SPECTRUM_PATH}/spectrum-mpi-rolling-release-gcc-4.9.3
export COMPILER_ROOT=${SPECTRUM_GCC}/bin
#export COMPILER_ROOT=/usr/tcetmp/bin

# Either of these pairs of BLAS/LAPACK libraries will work

export LAPACK_ROOT=${PKG_ROOT}/lapack/3.5.0
export BLAS_ROOT=${PKG_ROOT}/lapack/3.5.0
#export LAPACK_ROOT=/usr/tcetmp/packages/lapack/lapack-3.6.0-xlf-15.1.5
#export BLAS_ROOT=/usr/tcetmp/packages/blas/blas-3.6.0-xlf-15.1.5
export LAPACK_LDFLAGS="-L${LAPACK_ROOT}/lib -llapack -L${BLAS_ROOT}/lib -lblas"

./configure \
    --with-CXX=${COMPILER_ROOT}/mpicxx \
    --with-FC=/usr/tce/packages/gcc/gcc-4.9.3/bin/gfortran \
    --with-lapack=${LAPACK_ROOT} \
    --with-lapack-libs="${LAPACK_LDFLAGS}" \
    --with-hdf5=${PKG_ROOT}/hdf5/1.8.10p1 \
    --with-zlib=${PKG_ROOT}/zlib/1.2.3/lib \
    --with-gtest=no \
    --with-elemental=no \
    --enable-opt=yes \
    --enable-debug=no \
    "$@"

