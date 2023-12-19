# Toolchain for Lawrence Livermore National Laboratory CORAL machines
# (e.g., Lassen), assuming the following dependencies:
#
# - OpenMPI with GNU 12.2.1 compilers
# - LAPACK 3.11
# - HDF5 1.14.0

# Use MPI compiler wrappers because it simplifies detection of MPI
set(CMAKE_C_COMPILER /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-12.2.1/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-12.2.1/bin/mpicxx)
set(CMAKE_Fortran_COMPILER /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-12.2.1/bin/mpif90)

execute_process(
  COMMAND "echo $(which h5diff) | sed -e s/'/bin/h5diff'//"
  OUTPUT_VARIABLE hdf5_serial_path
)
set(HDF5_ROOT ${hdf5_serial_path})

set(BLAS_ROOT /usr/tcetmp/packages/lapack/lapack-3.11.0-gcc-11.2.1/)
set(LAPACK_ROOT /usr/tcetmp/packages/lapack/lapack-3.11.0-gcc-11.2.1/)

set(ScaLAPACK_ROOT /usr/tcetmp/packages/lapack/lapack-3.11.0-gcc-11.2.1/)
set(SCALAPACKDIR /usr/tcetmp/packages/lapack/lapack-3.11.0-gcc-11.2.1/)

set(CMAKE_EXE_LINKER_FLAGS_INIT "-lblas -llapack")

# Set flags to tune for the POWER9 architecture
set(CMAKE_C_FLAGS_INIT "-mcpu=powerpc64le -mtune=powerpc64le")
set(CMAKE_CXX_FLAGS_INIT "-mcpu=powerpc64le -mtune=powerpc64le")

set(CMAKE_LIBRARY_ARCHITECTURE ppc64le)
