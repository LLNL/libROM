#[[
This file is a part of ScalaWRAP.

Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
of the MIT license. See the top-level COPYRIGHT file for details.

SPDX-License-Identifier: MIT
]]
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules")

project("ScalaWRAP" LANGUAGES C CXX Fortran)

#[[-----------------------------------------------------------------------------
  | Set up the compiler and its flags. Besides the compiler the value the user
  | is most likely to want to give CMake is "OPTIM_FLAGS", which will replace
  | the release flags that CMake would use. For example, if building on the same
  | architecture where the code will run with the Intel compiler, a sensible
  | specification would be -DOPTIM_FLAGS="-O3 -xhost". Interprocedural
  | optimization is probably a good step to consider if you're serious about
  | performance.
  |
  | Right now only gcc, clang, or intel compilers are supported in this code.
  |--------------------------------------------------------------------------#]]

#[[-----------------------------------------------------------------------------
  | Detect the Fortran name mangling scheme and create a header with the macros
  | to do it.
  |--------------------------------------------------------------------------#]]
include(FortranCInterface)
FortranCInterface_HEADER(fc_interface.h MACRO_NAMESPACE "FC_")

#[[-----------------------------------------------------------------------------
  | Now detect ScaLAPACK. If MKL was found for LAPACK, by default we will use
  | the MKL ScaLAPACK. Otherwise, we look for -lscalapack. To override this
  | behavior, specify SCALAPACK_LIBNAMES. If it is not in the system's library
  | search path, you can specify locations to search in SCALAPACK_LIB_DIRS.
  |--------------------------------------------------------------------------#]]
if (NOT MKLROOT)
    set(MKLROOT $ENV{MKLROOT})
endif()
message(STATUS "MKLROOT = ${MKLROOT}")

if (NOT SCALAPACK_LIBNAMES)
    if (LAPACK_LIBRARIES MATCHES ".*mkl.*64.*" OR NOT ${MKLROOT} STREQUAL "")
        message(STATUS "Looking for MKL ScaLAPACK")
        find_library(SCALAPACK_LIBRARY NAMES mkl_scalapack_lp64)
        find_library(SCALAPACK_LIBRARY NAMES mkl_scalapack_lp64
            PATHS ${MKLROOT}/lib/intel64 ${SCALAPACK_LIB_DIRS})
        find_library(BLACS_LIBRARY NAMES mkl_blacs_intelmpi_lp64)
        find_library(BLACS_LIBRARY NAMES mkl_blacs_intelmpi_lp64
                     PATHS ${MKLROOT}/lib/intel64 ${SCALAPACK_LIB_DIRS})
    else()
        find_library(SCALAPACK_LIBRARY NAMES scalapack)
        find_library(SCALAPACK_LIBRARY NAMES scalapack PATHS ${SCALAPACK_LIB_DIRS})
    endif()
    set(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARY} ${BLACS_LIBRARY})
else()
    foreach(LIB ${SCALAPACK_LIBNAMES})
        find_library(SCALAPACK_LIBRARY NAMES ${LIB})
        find_library(SCALAPACK_LIBRARY NAMES ${LIB} PATHS ${SCALAPACK_LIB_DIRS})
        set(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARIES} ${SCALAPACK_LIBRARY})
    endforeach(LIB)
endif()
message(STATUS "ScaLAPACK libraries detected: ${SCALAPACK_LIBRARIES}")

add_library(scalawrap src/Context.cpp src/ScalaMat.cpp src/MPI_utils.cpp
    src/C_interface.cpp src/BLAS.cpp src/SVD.cpp)
target_include_directories(scalawrap
    PUBLIC
    ${MPI_C_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/src
    PRIVATE
    $<$<BOOL:${MKLROOT}>:${MKLROOT}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
target_link_libraries(scalawrap
    PUBLIC
    ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES} ${SCALAPACK_LIBRARIES}
    PRIVATE
)
if (SCALAPACK_LIBRARIES MATCHES ".*mkl.*")
    target_compile_definitions(scalawrap PRIVATE -DHAVE_MKL)
endif()

