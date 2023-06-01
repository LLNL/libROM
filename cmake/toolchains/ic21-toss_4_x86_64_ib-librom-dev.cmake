# Toolchain for Lawrence Livermore National Laboratory TOSS4 machines
# (e.g., Quartz), assuming the following dependencies:
#
# - MVAPICH2/2.3 with Intel 21.6.0 compilers
# - Intel MKL
# - HDF5 1.14.0

# Use MPI compiler wrappers because it simplifies detection of MPI
set(CMAKE_C_COMPILER /usr/tce/packages/mvapich2/mvapich2-2.3.7-intel-classic-2021.6.0-magic/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/tce/packages/mvapich2/mvapich2-2.3.7-intel-classic-2021.6.0-magic/bin/mpicxx)
set(CMAKE_Fortran_COMPILER /usr/tce/packages/mvapich2/mvapich2-2.3.7-intel-classic-2021.6.0-magic/bin/mpif90)
set(BLA_VENDOR Intel10_64lp)
set(HDF5_ROOT /usr/tce/packages/hdf5/hdf5-1.14.0-mvapich2-2.3.7-intel-classic-2021.6.0/)

set(BLAS_ROOT /usr/tce/packages/mkl/mkl-2022.1.0)
set(LAPACK_ROOT /usr/tce/packages/mkl/mkl-2022.1.0)
set(CMAKE_LIBRARY_ARCHITECTURE intel64)

set(MKL_MESSAGE_PASSING MPICH2)
set(MKL_THREADING Sequential)
set(MKL_ROOT /usr/tce/packages/mkl/mkl-2022.1.0)
