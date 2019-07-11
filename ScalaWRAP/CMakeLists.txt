#[[
This file is a part of ScalaWRAP.

Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
of the MIT license. See the top-level COPYRIGHT file for details.

SPDX-License-Identifier: MIT
]]

#[[-----------------------------------------------------------------------------
  | General CMake setup stuff - set language standards, CMake version, build
  | type, project info.
  |--------------------------------------------------------------------------#]]

cmake_minimum_required(VERSION 3.14.5)

# C++11 is the minimum standard version supported.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set (CMAKE_C_STANDARD 99)

# I do not want non-standard code to work in this library.
set(CMAKE_CXX_EXTENSIONS OFF)

# Build an optimized library by default.
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No CMAKE_BUILD_TYPE specified; configuring for Release")
    set(CMAKE_BUILD_TYPE Release)
endif()

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

if (OPTIM_FLAGS)
    message(STATUS "Using user specified optimization flags for C/C++: ${OPTIM_FLAGS}")
endif()

# These options are compatible with gcc, clang, and Intel.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_DEBUG "-g")

if (OPTIM_FLAGS)
    set(CMAKE_C_FLAGS_RELEASE ${OPTIM_FLAGS})
    set(CMAKE_CXX_FLAGS_RELEASE ${OPTIM_FLAGS})
else()
    set(CMAKE_C_FLAGS_RELEASE "-O3")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

#[[-----------------------------------------------------------------------------
  | Detect the Fortran name mangling scheme and create a header with the macros
  | to do it.
  |--------------------------------------------------------------------------#]]
include(FortranCInterface)
FortranCInterface_HEADER(fc_interface.h MACRO_NAMESPACE "FC_")

#[[-----------------------------------------------------------------------------
  | Find the MPI vendor automatically. This should set everything up correctly
  |--------------------------------------------------------------------------#]]
find_package(MPI REQUIRED)

#[[-----------------------------------------------------------------------------
  | Find LAPACK. I start by looking for the MKL. This library assumes a
  | 32-bit integer interface so we only look for that one. Variables to control
  | this detection: set LAPACK_VENDOR to a specific vendor. See CMake docs on
  | FindLAPACK for valid vendors but note that it should be be the lp64 model
  | and not the ilp64 model if linking against the MKL.
  |
  | Set LAPACK_EXECUTION to one of "sequential" or "threaded" to override that
  | preference.
  |
  | By default, I look for sequential MKL first, then threaded, then search all
  | libraries. If the detected LAPACK is the ilp64 interface to MKL, this is a
  | fatal error. Sequential is preferred by default because it is assumed that
  | you will run a process per core for parallelism; if you would rather have
  | a multithreaded model set LAPACK_EXECUTION to "threaded"
  |--------------------------------------------------------------------------#]]

# If someone's read the CMake docs they might try to set the vendor with
# BLA_VENDOR; I added LAPACK_VENDOR as the variable because really who shortens
# it to BLA?
if (BLA_VENDOR AND (NOT LAPACK_VENDOR))
    set(LAPACK_VENDOR ${BLA_VENDOR})
endif()
if (LAPACK_VENDOR)
    set(BLA_VENDOR ${LAPACK_VENDOR})
    message(STATUS "Looking for LAPACK w/ vendor = ${BLA_VENDOR}")
    find_package(LAPACK REQUIRED)
elseif(LAPACK_EXECUTION)
    if (LAPACK_EXECUTION STREQUAL threaded)
        set(BLA_VENDOR Intel10_64lp)
        message(STATUS "Looking for LAPACK w/ vendor = ${BLA_VENDOR}")
        find_package(LAPACK)
        if (NOT LAPACK_FOUND)
            unset(BLA_VENDOR)
            message(STATUS "Looking for LAPACK from any vendor")
            find_package(LAPACK REQUIRED)
        endif()
    else()
        set(BLA_VENDOR Intel10_64lp_seq)
        message(STATUS "Looking for LAPACK w/ vendor = ${BLA_VENDOR}")
        find_package(LAPACK )
        if (NOT LAPACK_FOUND)
            unset(BLA_VENDOR)
            message(STATUS "Looking for LAPACK from any vendor")
            find_package(LAPACK REQUIRED)
        endif()
    endif()
else()
    set(BLA_VENDOR Intel10_64lp_seq)
    message(STATUS "Looking for LAPACK w/ vendor = ${BLA_VENDOR}")
    find_package(LAPACK)
    if (NOT LAPACK_FOUND)
        set(BLA_VENDOR Intel10_64lp)
        message(STATUS "Looking for LAPACK w/ vendor = ${BLA_VENDOR}")
        find_package(LAPACK)
        if (NOT LAPACK_FOUND)
            set(BLA_VENDOR OpenBLAS)
            message(STATUS "Looking for LAPACK w/ vendor = ${BLA_VENDOR}")
            find_package(LAPACK)
            if (NOT LAPACK_FOUND)
                unset(BLA_VENDOR)
                message(STATUS "Looking for LAPACK from any vendor")
                find_package(LAPACK REQUIRED)
            endif()
        endif()
    endif()
endif()

if (LAPACK_LIBRARIES MATCHES ".*ilp64.*")
    message(FATAL_ERROR "Cannot link against ILP64 interface to MKL")
endif()

#[[-----------------------------------------------------------------------------
  | Now detect ScaLAPACK. If MKL was found for LAPACK, by default we will use
  | the MKL ScaLAPACK. Otherwise, we look for -lscalapack. To override this
  | behavior, specify SCALAPACK_LIBNAMES. If it is not in the system's library
  | search path, you can specify locations to search in SCALAPACK_LIB_DIRS.
  |--------------------------------------------------------------------------#]]
if (NOT MKLROOT)
    set(MKLROOT $ENV{MKLROOT})
endif()

if (NOT SCALAPACK_LIBNAMES)
    if (LAPACK_LIBRARIES MATCHES ".*mkl.*64.*")
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

add_subdirectory(test)
