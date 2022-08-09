###############################################################################
#
#  Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
#  and other libROM project developers. See the top-level COPYRIGHT
#  file for details.
#
#  SPDX-License-Identifier: (Apache-2.0 OR MIT)
#
###############################################################################

# Require CMake version 3.12 or later to expose modern features of
# CMake for maintainability. With additional work, it is possible to
# rewrite the build system generator for CMake 3.10+, 3.8+, 3.1+,
# 3.0+, or 2.8+, with increasing amounts of work as the minimum
# required version is decreased. Notably, CMake's find modules export
# targets (e.g., FindZLIB.cmake exports the ZLIB::ZLIB IMPORTED
# target) as of version 3.1, and CMake supports using the
# <PackageName>_ROOT variables to set the prefix search path for
# find_package(<PackageName>) as of version 3.12.
cmake_minimum_required(VERSION 3.12)
cmake_policy(SET CMP0074 NEW) # Use <PackageName>_ROOT variables in find_package

# Even though this project is a pure C++ library, the C language must be enabled
# for MPI and HDF5 dependency checking. The Fortran language must be enabled
# for Fortran-C name mangling convention detection.
project(libROM
  VERSION 1.0.0
  DESCRIPTION "Model reduction library emphasizing large-scale parallelism"
  HOMEPAGE_URL "https://github.com/LLNL/libROM"
  LANGUAGES C CXX Fortran)

set(libROM_CMAKE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(libROM_CMAKE_MODULE_PATH ${libROM_CMAKE_PATH}/modules)
list(APPEND CMAKE_MODULE_PATH ${libROM_CMAKE_MODULE_PATH})

option(USE_MFEM "Build libROM with MFEM" OFF)
option(BUILD_STATIC "Build libROM as a static library" OFF)

## Set a bunch of variables to generate a configure header
# Enable assertion checking if debug symbols generated
if((CMAKE_BUILD_TYPE STREQUAL "Debug") OR
    (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo"))
  set(DEBUG_CHECK_ASSERTIONS "1")
endif((CMAKE_BUILD_TYPE STREQUAL "Debug") OR
    (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo"))

set(CAROM_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(CAROM_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(CAROM_VERSION_PATCHLEVEL "${PROJECT_VERSION_PATCH}")
set(CAROM_PACKAGE ${CMAKE_PROJECT_NAME})
set(CAROM_PACKAGE_NAME ${CMAKE_PROJECT_NAME})
set(CAROM_PACKAGE_STRING "${CMAKE_PROJECT_NAME} ${CMAKE_PROJECT_VERSION}")
set(CAROM_PACKAGE_TARNAME ${CMAKE_PROJECT_NAME})
set(CAROM_PACKAGE_URL ${CMAKE_PROJECT_HOMEPAGE_URL})
set(CAROM_PACKAGE_VERSION ${CMAKE_PROJECT_VERSION})
set(CAROM_VERSION ${CMAKE_PROJECT_VERSION})

include(CheckIncludeFiles)
check_include_files(dlfcn.h CAROM_HAVE_DLFCN_H)
check_include_files(inttypes.h CAROM_HAVE_INTTYPES_H)
check_include_files(memory.h CAROM_HAVE_MEMORY_H)
check_include_files(stdint.h CAROM_HAVE_STDINT_H)
check_include_files(stdlib.h CAROM_HAVE_STDLIB_H)
check_include_files(strings.h CAROM_HAVE_STRINGS_H)
check_include_files(string.h CAROM_HAVE_STRING_H)
check_include_files(sys/stat.h CAROM_HAVE_SYS_STAT_H)
check_include_files(sys/types.h CAROM_HAVE_SYS_TYPES_H)
check_include_files(unistd.h CAROM_HAVE_UNISTD_H)

# Construct an expression consisting of the 24 ANSI C headers
set(stdc_header_list "assert.h;complex.h;ctype.h;errno.h;fenv.h;float.h")
set(stdc_header_list "${stdc_header_list};inttypes.h;iso646.h;limits.h")
set(stdc_header_list "${stdc_header_list};locale.h;math.h;setjmp.h;signal.h")
set(stdc_header_list "${stdc_header_list};stdarg.h;stdbool.h;stdint.h")
set(stdc_header_list "${stdc_header_list};stddef.h;stdio.h;stdlib.h;string.h")
set(stdc_header_list "${stdc_header_list};tgmath.h;time.h;wchar.h;wctype.h")
check_include_files("${stdc_header_list}" CAROM_STDC_HEADERS)

# Define variables for use in generating a configure file
if(GTEST_FOUND)
  set(CAROM_HAS_GTEST 1)
endif(GTEST_FOUND)

if(BLAS_FOUND)
  set(CAROM_HAVE_BLAS 1)
endif(BLAS_FOUND)

if(LAPACK_FOUND)
  set(CAROM_HAVE_LAPACK 1)
endif(LAPACK_FOUND)

if(HDF5_FOUND)
  set(CAROM_HAVE_HDF5 1)
endif(HDF5_FOUND)

# List minimum version requirements for dependencies where possible to make
# packaging easier later.
find_package(HDF5 1.8.0 REQUIRED)

find_package(BLAS 3.4.0 REQUIRED)
find_package(LAPACK 3.4.0 REQUIRED)

# In FindMPI.cmake, "CXX" refers to "the MPI C API [being] usable from C++"
find_package(MPI 1.2 REQUIRED)

find_package(ZLIB 1.2.3 REQUIRED)

find_package(Doxygen 1.8.5)

find_package(GTest 1.6.0)

if (USE_MFEM)
    find_library(MFEM mfem "${CMAKE_SOURCE_DIR}/dependencies/mfem")
    find_library(HYPRE HYPRE "${CMAKE_SOURCE_DIR}/dependencies/hypre/src/hypre/lib")
    find_library(PARMETIS parmetis "${CMAKE_SOURCE_DIR}/dependencies/parmetis-4.0.3/build/lib/libparmetis")
    find_library(METIS metis "${CMAKE_SOURCE_DIR}/dependencies/parmetis-4.0.3/build/lib/libmetis")
    find_path(MFEM_INCLUDES mfem.hpp "${CMAKE_SOURCE_DIR}/dependencies/mfem")
    find_path(HYPRE_INCLUDES HYPRE.h "${CMAKE_SOURCE_DIR}/dependencies/hypre/src/hypre/include")
    find_path(PARMETIS_INCLUDES metis.h "${CMAKE_SOURCE_DIR}/dependencies/parmetis-4.0.3/metis/include")
endif()

add_subdirectory(lib)

# Use the C++11 standard as an entire feature instead of
# enumerating individual compiler features for simplicity
target_compile_features(ROM PRIVATE cxx_std_11)

if (USE_MFEM)
  set(examples
    poisson_global_rom
    poisson_local_rom_greedy
    dg_advection_global_rom
    dg_advection_local_rom_matrix_interp
    mixed_nonlinear_diffusion
    nonlinear_elasticity_global_rom
    dg_advection
    nonlinear_elasticity
    heat_conduction
    parametric_heat_conduction
    dg_euler
    local_tw_csv
    parametric_tw_csv)
  set(example_directories
    prom
    prom
    prom
    prom
    prom
    prom
    dmd
    dmd
    dmd
    dmd
    dmd
    dmd
    dmd)

  list(LENGTH examples len1)
  math(EXPR len2 "${len1} - 1")

  foreach(val RANGE ${len2})
    list(GET examples ${val} name)
    list(GET example_directories ${val} example_dir)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/${example_dir})
    add_executable(${name} examples/${example_dir}/${name}.cpp)

    target_link_libraries(${name}
      PRIVATE ROM ${MPI_C_LINK_FLAGS} ${MPI_C_LIBRARIES} MPI::MPI_C ${MPI_FORTRAN_LINK_FLAGS} ${MPI_FORTRAN_LIBRARIES} MPI::MPI_Fortran)
    target_include_directories(${name}
      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
      ${MPI_C_INCLUDE_DIRS})
    target_compile_features(${name} PRIVATE cxx_std_11)
  endforeach() # IN LISTS examples
  file(COPY examples/data DESTINATION ${CMAKE_BINARY_DIR}/examples)
  file(COPY examples/dmd/heat_conduction_csv.sh DESTINATION ${CMAKE_BINARY_DIR}/examples/dmd)
endif()

set(misc_example_names
  combine_samples)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/misc)

foreach(name IN LISTS misc_example_names)
  add_executable(${name} examples/misc/${name}.cpp)

  target_link_libraries(${name}
    PRIVATE ROM ${MPI_C_LINK_FLAGS} ${MPI_C_LIBRARIES} MPI::MPI_C ${MPI_FORTRAN_LINK_FLAGS} ${MPI_FORTRAN_LIBRARIES} MPI::MPI_Fortran)
  target_include_directories(${name}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${MPI_C_INCLUDE_DIRS})
  target_compile_features(${name} PRIVATE cxx_std_11)
endforeach(name) # IN LISTS misc_exmaple_names

set(regression_test_names
  smoke_test
  test_include
  uneven_dist
  weak_scaling
  random_test
  smoke_static
  load_samples)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

foreach(name IN LISTS regression_test_names)
  add_executable(${name} tests/${name}.cpp)

  target_link_libraries(${name}
    PRIVATE ROM ${MPI_C_LINK_FLAGS} ${MPI_C_LIBRARIES} MPI::MPI_C ${MPI_FORTRAN_LINK_FLAGS} ${MPI_FORTRAN_LIBRARIES} MPI::MPI_Fortran)
  target_include_directories(${name}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${MPI_C_INCLUDE_DIRS})
  target_compile_features(${name} PRIVATE cxx_std_11)
endforeach(name) # IN LISTS regression_test_names

if(GTEST_FOUND)
  set(unit_test_stems
    Vector
    Matrix
    DEIM
    DMD
    GNAT
    QDEIM
    S_OPT
    SVD
    StaticSVD
    RandomizedSVD
    IncrementalSVD
    GreedyCustomSampler)
  foreach(stem IN LISTS unit_test_stems)
    add_executable(test_${stem} tests/test_${stem}.cpp)
    target_link_libraries(test_${stem} PRIVATE ROM
      ${MPI_C_LINK_FLAGS} ${MPI_C_LIBRARIES} MPI::MPI_C ${MPI_FORTRAN_LINK_FLAGS} ${MPI_FORTRAN_LIBRARIES} MPI::MPI_Fortran GTest::GTest)
    target_include_directories(test_${stem}
      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
      ${MPI_C_INCLUDE_DIRS})
    target_compile_features(test_${stem} PRIVATE cxx_std_11)
    target_compile_definitions(test_${stem} PRIVATE CAROM_HAS_GTEST)
  endforeach(stem) # IN LISTS unit_test_stems
endif(GTEST_FOUND)

# NOTE(goxberry@gmail.com, oxberry1@llnl.gov): This code snippet
# builds the Doxygen documentation, but outputs said documentation to
# ${CMAKE_CURRENT_SOURCE_DIR}/docs)...
if(DOXYGEN_FOUND)
  set(doxyfile ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile)

  add_custom_target(
    documentation ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM
    )
  add_dependencies(documentation ROM)

  add_custom_target(
    doxygen_tagfile
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/rom.tag
    ${CMAKE_CURRENT_BINARY_DIR}/docs/html/rom.tag)
  add_dependencies(doxygen_tagfile documentation)

endif(DOXYGEN_FOUND)
