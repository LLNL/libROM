# Toolchain assuming a macOS system with all prerequisites installed
# via Homebrew (https://brew.sh/).
#
# brew install openmpi && brew unlink openmpi
# brew install mpich && brew unlink mpich && brew link openmpi
# brew install openblas
# brew install hdf5
# brew install szip

# Use MPI compiler wrappers because it simplifies detection of MPI
set(CMAKE_C_COMPILER /usr/local/opt/mpich/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/local/opt/mpich/bin/mpicxx)
set(CMAKE_Fortran_COMPILER /usr/local/opt/mpich/bin/mpif90)
set(MPI_ROOT /usr/local/opt/mpich)

set(BLA_VENDOR Intel10_64lp_seq)
set(BLAS_ROOT /opt/intel/mkl)
set(LAPACK_ROOT /opt/intel/mkl)
set(MKL_ROOT /opt/intel/mkl)
set(MKL_MESSAGE_PASSING MPICH)
set(MKL_THREADING Sequential)

set(HDF5_ROOT /usr/local/opt/hdf5)
#string(APPEND CMAKE_CXX_FLAGS_INIT " -gxx-name=/usr/tce/packages/gcc/gcc-4.9.3/bin")
