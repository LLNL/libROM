###############################################################################
#
#  Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
#  and other libROM project developers. See the top-level COPYRIGHT
#  file for details.
#
#  SPDX-License-Identifier: (Apache-2.0 OR MIT)
#
###############################################################################

# Debugging tip: If using the GNU Make generator (the default), build
# with `make VERBOSE=1` to see flags being used by GNU Make. This
# output helps in trying to diagnose configuration errors. The
# equivalent for the Ninja generator is `ninja -v`. Also note that
# because this project "uses Fortran" -- a Fortran compiler is
# required to infer Fortran name-mangling conventions -- KitWare's
# "ninja-fortran" is needed to run the resulting Ninja files. This
# software can be installed, e.g., via Spack (`spack install
# ninja-fortran`).

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
  VERSION 0.0.1
  DESCRIPTION "Model reduction library emphasizing large-scale parallelism"
  HOMEPAGE_URL "https://github.com/LLNL/libROM"
  LANGUAGES C CXX Fortran)

set(libROM_CMAKE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(libROM_CMAKE_MODULE_PATH ${libROM_CMAKE_PATH}/modules)
list(APPEND CMAKE_MODULE_PATH ${libROM_CMAKE_MODULE_PATH})

option(USE_MFEM "Build libROM with MFEM" OFF)
option(BUILD_STATIC "Build libROM as a static library" OFF)

include(FortranCInterface)
FortranCInterface_HEADER(${CMAKE_SOURCE_DIR}/FCMangle.h
  MACRO_NAMESPACE "CAROM_FC_")
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

configure_file(${CMAKE_SOURCE_DIR}/CAROM_config.h.in
  ${CMAKE_SOURCE_DIR}/CAROM_config.h @ONLY)

# While it is tempting to use file globbing here, file globbing is
# considered a "modern CMake anti-pattern" because "CMake is not a
# build system -- it is a build system _generator_". Instead, use
# some file name "stems" for generating file names. This construction
# is useful when files may be moved to different directories.
set(module_list
  Matrix
  Vector
  BasisGenerator
  BasisReader
  BasisWriter
  svd/IncrementalSVD
  svd/IncrementalSVDFastUpdate
  svd/IncrementalSVDStandard
  svd/RandomizedSVD
  svd/StaticSVD
  svd/SVD
  greedy/GreedyCustomSampler
  greedy/GreedyRandomSampler
  greedy/GreedySampler
  manifold_interp/Interpolater
  manifold_interp/MatrixInterpolater
  manifold_interp/VectorInterpolater
  DMD
  hyperreduction/DEIM
  hyperreduction/GNAT
  hyperreduction/QDEIM
  hyperreduction/STSampling
  utils/Database
  utils/HDFDatabase
  utils/Utilities
  utils/ParallelBuffer)
set(source_files)
foreach(module IN LISTS module_list)
  list(APPEND source_files ${module}.C)
endforeach(module) # IN LISTS module_list

list(APPEND source_files
  Options.h)

if (USE_MFEM)
    list(APPEND source_files
      mfem/SampleMesh.hpp)
endif()

list(APPEND source_files
  scalapack_c_wrapper.c
  scalapack_f_wrapper.f90)
if (BUILD_STATIC)
  add_library(ROM ${source_files})
else()
  add_library(ROM SHARED ${source_files})
endif()

# List minimum version requirements for dependencies where possible to make
# packaging easier later.
find_package(HDF5 1.8.0 REQUIRED)

find_package(BLAS 3.4.0 REQUIRED)
find_package(LAPACK 3.4.0 REQUIRED)

# If MKL libraries not found, search for reference ScaLAPACK. If MKL
# libraries found, search for MKL ScaLAPACK; if MKL ScaLAPACK not
# found, search for reference ScaLAPACK. It seems that only
if (BLAS_LIBRARIES MATCHES ".*mkl.*")
  find_package(MKL COMPONENTS BLACS ScaLAPACK)
  if (NOT MKL_ScaLAPACK_FOUND)
    find_package(ScaLAPACK REQUIRED)
    target_link_libraries(ROM PUBLIC ${ScaLAPACK_LIBRARIES})
  else()
    target_link_libraries(ROM PUBLIC ${MKL_ScaLAPACK_LIBRARY} ${MKL_BLACS_LIBRARY} ${MKL_LIBRARIES})
    target_include_directories(ROM PUBLIC ${MKL_INCLUDE_DIRS})
  endif()
else() # BLAS or LAPACK isn't MKL
  find_package(ScaLAPACK REQUIRED)
  target_link_libraries(ROM PUBLIC ${ScaLAPACK_LIBRARIES})
endif()

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

# PUBLIC dependencies are transitive; these dependencies are used in
# the API headers *and* in their implementations
#
# INTERFACE dependencies are used in the API headers, but not in the
# API implementation (e.g., API forwards objects to another library
# without dereferencing those objects -- only pointers are used)
#
# PRIVATE dependencies are used in the API implementation, but not in
# the headers
#
# Using both the MPI::MPI_C target and the two MPI_C "flag variables
# (i.e., "MPI_C_LINK_FLAGS, MPI_C_LIBRARIES) is probably redundant,
# but is done here to ease a potential rollback to CMake 2.8 or CMake
# 3.0.
target_link_libraries(ROM
  PUBLIC ${MPI_C_LINK_FLAGS} ${MPI_C_LIBRARIES} MPI::MPI_C ${MPI_FORTRAN_LINK_FLAGS} ${MPI_FORTRAN_LIBRARIES} MPI::MPI_Fortran ${HDF5_LIBRARIES}
  ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${MFEM} ${HYPRE} ${PARMETIS} ${METIS}
  PRIVATE ${ZLIB_LIBRARIES} ZLIB::ZLIB)

target_include_directories(ROM PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${MFEM_INCLUDES}
  ${HYPRE_INCLUDES}
  ${PARMETIS_INCLUDES}
  ${HDF5_C_INCLUDE_DIRS}
  ${MPI_C_INCLUDE_DIRS}
  ${MFEM_C_INCLUDE_DIRS}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Use the C++11 standard as an entire feature instead of
# enumerating individual compiler features for simplicity
target_compile_features(ROM PRIVATE cxx_std_11)

if (USE_MFEM)
  set(examples
    poisson_global_rom
    poisson_local_rom_greedy
    dg_advection_global_rom
    dg_advection_local_rom_matrix_interp
    dg_advection
    nonlinear_elasticity
    heat_conduction
    dg_euler)
  set(example_directories
    prom
    prom
    prom
    prom
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
endif()

set(regression_test_names
  smoke_test
  uneven_dist
  weak_scaling
  random_test
  smoke_static
  load_samples)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

foreach(name IN LISTS regression_test_names)
  add_executable(${name} tests/${name}.C)

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
    SVD
    StaticSVD
    RandomizedSVD
    IncrementalSVD
    GreedyCustomSampler)
  foreach(stem IN LISTS unit_test_stems)
    add_executable(test_${stem} tests/test_${stem}.C)
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
