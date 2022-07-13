# Toolchain for Lawrence Livermore National Laboratory TOSS3 machines
# (e.g., Quartz), assuming the following dependencies:
#
# - MVAPICH2/2.2 with GNU 4.9.3 compilers
# - Intel MKL
# - HDF5 1.8.17

# Use MPI compiler wrappers because it simplifies detection of MPI
set(CMAKE_C_COMPILER /usr/tce/packages/mvapich2/mvapich2-2.2-gcc-4.9.3/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/tce/packages/mvapich2/mvapich2-2.2-gcc-4.9.3/bin/mpicxx)
set(CMAKE_Fortran_COMPILER /usr/tce/packages/mvapich2/mvapich2-2.2-gcc-4.9.3/bin/mpif90)
set(BLA_VENDOR Intel10_64lp)
set(HDF5_ROOT /usr/tce/packages/hdf5/hdf5-serial-1.8.17-gcc4.9.3)
#string(APPEND CMAKE_CXX_FLAGS_INIT " -gxx-name=/usr/tce/packages/gcc/gcc-7.1.0/bin/g++")

set(BLAS_ROOT /usr/tce/packages/mkl/mkl-2019.0/mkl)
set(LAPACK_ROOT /usr/tce/packages/mkl/mkl-2019.0/mkl)
set(CMAKE_LIBRARY_ARCHITECTURE intel64)

set(MKL_MESSAGE_PASSING MPICH2)
set(MKL_THREADING Sequential)
set(MKL_ROOT /usr/tce/packages/mkl/mkl-2019.0/mkl)
