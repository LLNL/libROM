# Toolchain for Lawrence Livermore National Laboratory TOSS3 machines
# (e.g., Quartz), using the compilers and HDF5 library as available
# in the environment (i.e., through "module load").
#
# Note that:
#   which mpicc
#   which mpicxx
#   which mpif90
#   which h5diff
# should point to MPI compilers and the HDF5 library built with
# the same C/FORTRAN compilers.

# Use MPI compiler wrappers because it simplifies detection of MPI
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_Fortran_COMPILER mpif90)
set(CMAKE_Fortran_FLAGS -mcmodel=large)
#set(BLA_VENDOR Intel10_64lp)

execute_process(
  COMMAND "echo $(which h5diff) | sed -e s/'/bin/h5diff'//"
  OUTPUT_VARIABLE hdf5_serial_path
)
set(HDF5_ROOT ${hdf5_serial_path})

#set(BLAS_ROOT /usr/tce/packages/mkl/mkl-2019.0/mkl)
#set(LAPACK_ROOT /usr/tce/packages/mkl/mkl-2019.0/mkl)
set(BLAS_ROOT /usr/lib/x86_64-linux-gnu/blas)
set(LAPACK_ROOT /usr/lib/x86_64-linux-gnu/blas)

set(CMAKE_LIBRARY_ARCHITECTURE intel64)

set(MKL_MESSAGE_PASSING MPICH2)
set(MKL_THREADING Sequential)
set(MKL_ROOT /usr/tce/packages/mkl/mkl-2019.0/mkl)
