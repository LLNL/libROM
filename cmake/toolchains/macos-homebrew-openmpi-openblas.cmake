# Toolchain assuming a macOS system with all prerequisites installed
# via Homebrew (https://brew.sh/).
#
# brew install openmpi
# brew install openblas
# brew install hdf5
# brew install szip
# brew install scalapack

# Use MPI compiler wrappers because it simplifies detection of MPI
set(MPI_ROOT /usr/local/opt/openmpi)
set(CMAKE_C_COMPILER /usr/local/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/local/bin/mpicxx)
set(CMAKE_Fortran_COMPILER /usr/local/bin/mpif90)

set(BLA_VENDOR OpenBLAS)
set(BLAS_ROOT /usr/local/opt/openblas)
set(LAPACK_ROOT /usr/local/opt/openblas)

set(ScaLAPACK_ROOT /usr/local/opt/scalapack)

set(HDF5_ROOT /usr/local/opt/hdf5)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Output compilation database" FORCE)
