# Prefix directories for various library installations
execute_process(COMMAND spack location -i hdf5 OUTPUT_VARIABLE
  HDF5_ROOT OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND spack location -i openblas OUTPUT_VARIABLE
  BLAS_ROOT OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND spack location -i openblas OUTPUT_VARIABLE
  LAPACK_ROOT OUTPUT_STRIP_TRAILING_WHITESPACE)
set(BLA_VENDOR OpenBLAS)

execute_process(COMMAND spack location -i openmpi OUTPUT_VARIABLE
  MPI_ROOT OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CMAKE_CXX_COMPILER "${MPI_ROOT}/bin/mpicxx")
set(CMAKE_C_COMPILER "${MPI_ROOT}/bin/mpicc")
set(CMAKE_Fortran_COMPILER "${MPI_ROOT}/bin/mpif90")

execute_process(COMMAND spack location -i netlib-scalapack^openmpi OUTPUT_VARIABLE
  ScaLAPACK_ROOT OUTPUT_STRIP_TRAILING_WHITESPACE)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Output compilation database" FORCE)
