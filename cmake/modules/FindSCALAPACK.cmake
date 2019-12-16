# CMake "find module" for ScaLAPACK. Even though reference ScaLAPACK
# provides a CMake config file, the Intel MKL version of ScaLAPACK
# does not. This CMake find module can be bypassed in favor of
# detecting the reference ScaLAPACK find module via setting
# CMAKE_FIND_PACKAGE_PREFER_CONFIG to TRUE.

# General strategy:
#
# 1) Easiest: check to see if reference ScaLAPACK installed.
#
# 2) Check to see if PBLAS is needed. (If it is, fail for now, then
# implement a FindPBLAS.cmake later on.)
#
# 3) Check to see if BLACS is needed. (If it is, fail for now, then
# implement a FindBLACS.cmake later on.)
#
# 4) Probably need to implement the basics of a FindMKL.cmake package.
#

# Search hardcoded paths as follows:
# - /usr/lib64 and /usr/lib for system-wide Linux/BSD installation via system
#   package manager
# - /usr/local/lib64 and /usr/local/lib for system-wide Linux/BSD installation
#    outside of system package manager; includes Homebrew
# - /opt/local for MacPorts installation
# - /sw/lib or /opt/sw/lib for Fink installation
find_library(ScaLAPACK_LIBRARY
  NAMES scalapack scalapack-mpi scalapack-mpich scalapack-mpich2 scalapack-openmpi
  PATHS /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib
  /opt/local/lib /opt/sw/lib /sw/lib
  ENV LD_LIBRARY_PATH
  ENV DYLD_FALLBACK_LIBRARY_PATH
  ENV DYLD_LIBRARY_PATH
  ENV SCALAPACKDIR
  ENV BLACSDIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ScaLAPACK
  REQUIRED_VARS ScaLAPACK_LIBRARY
  )
set(ScaLAPACK_LIBRARIES ${ScaLAPACK_LIBRARY})
