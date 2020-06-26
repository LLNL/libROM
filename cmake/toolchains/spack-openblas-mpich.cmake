# Prefix directories for various library installations
execute_process(COMMAND spack location -i hdf5 OUTPUT_VARIABLE
  HDF5_ROOT OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND spack location -i openblas OUTPUT_VARIABLE
  BLAS_ROOT OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND spack location -i openblas OUTPUT_VARIABLE
  LAPACK_ROOT OUTPUT_STRIP_TRAILING_WHITESPACE)
set(BLA_VENDOR OpenBLAS)

execute_process(COMMAND spack location -i mpich OUTPUT_VARIABLE
  MPI_ROOT OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CMAKE_CXX_COMPILER "${MPI_ROOT}/bin/mpicxx")
set(CMAKE_C_COMPILER "${MPI_ROOT}/bin/mpicc")
set(CMAKE_Fortran_COMPILER "${MPI_ROOT}/bin/mpif90")

# NOTE(goxberry@gmail.com, oxberry1@llnl.gov): The version specifier
# 2.0.2 is necessary as of 2019-12-17 because the ^mpich flag in the
# netlib-scalapack spack spec does not by itself uniquely identify the
# ScaLAPACK installation on a machine that uses macOS Homebrew,
# installs netlib-scalapack and open-mpi via Homebrew, registers these
# installations in spack, then installs mpich and
# netlib-scalapack^mpich via spack. Strangely, `spack install
# netlib-scalapack^mpich` uniquely identifies a spack concretization,
# but `spack location -i netlib-scalapack^mpich` does not.

execute_process(COMMAND spack location -i netlib-scalapack@2.0.2^mpich OUTPUT_VARIABLE
  ScaLAPACK_ROOT OUTPUT_STRIP_TRAILING_WHITESPACE)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Output compilation database" FORCE)
