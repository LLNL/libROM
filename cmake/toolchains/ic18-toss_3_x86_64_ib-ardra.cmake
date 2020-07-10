# Toolchain for Lawrence Livermore National Laboratory TOSS3 machines
# (e.g., Quartz), assuming the following dependencies:
#
# - MVAPICH2 with Intel 18.0.2 compilers
# - Intel MKL
# - HDF5 1.8.10p1

# Use MPI compiler wrappers because it simplifies detection of MPI
set(CMAKE_C_COMPILER /usr/tce/packages/mvapich2/mvapich2-2.2-intel-18.0.2/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/tce/packages/mvapich2/mvapich2-2.2-intel-18.0.2/bin/mpicxx)
set(CMAKE_Fortran_COMPILER /usr/tce/packages/mvapich2/mvapich2-2.2-intel-18.0.2/bin/mpif90)
set(BLA_VENDOR Intel10_64lp)
set(HDF5_ROOT /usr/gapps/bdiv/toss_3_x86_64_ib/intel-18-mvapich2-2.2/hdf5/1.8.10p1)
string(APPEND CMAKE_CXX_FLAGS_INIT " -gxx-name=/usr/tce/packages/gcc/gcc-7.1.0/bin/g++")

set(BLAS_ROOT /usr/tce/packages/mkl/mkl-2019.0/mkl)
set(LAPACK_ROOT /usr/tce/packages/mkl/mkl-2019.0/mkl)
set(CMAKE_LIBRARY_ARCHITECTURE intel64)

set(MKL_MESSAGE_PASSING MPICH2)
set(MKL_THREADING Sequential)
set(MKL_ROOT /usr/tce/packages/mkl/mkl-2019.0/mkl)
