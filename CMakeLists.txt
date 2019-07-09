# Require CMake version 3.10 or later to expose modern features of
# CMake for maintainability. With additional work, it is possible to
# rewrite the build system generator for CMake 3.8+, 3.1+, 3.0+, or
# 2.8+, with increasing amounts of work as the minimum required
# version is decreased.
cmake_minimum_required(VERSION 3.10)

project(libROM
  VERSION 0.0.1
  DESCRIPTION "Model reduction library emphasizing large-scale parallelism"
  HOMEPAGE_URL "https://github.com/LLNL/libROM"
  LANGUAGES C CXX)

# List minimum version requirements for dependencies where possible to make
# packaging easier later.
find_package(HDF5 1.8.0 REQUIRED)
find_package(BLAS 3.4.0 REQUIRED)
find_package(LAPACK 3.4.0 REQUIRED)

# In FindMPI.cmake, "CXX" refers to "the MPI C API [being] usable from C++"
find_package(MPI 1.2 REQUIRED)

find_package(ZLIB 1.2.3 REQUIRED)

find_package(Doxygen 1.8.10)

# find_package(Perl 5.16.3 REQUIRED)

find_package(GTest 1.6.0)

# While it is tempting to use file globbing here, file globbing is
# considered a "modern CMake anti-pattern" because "CMake is not a
# build system -- it is a build system _generator_". Instead, use
# some file name "stems" for generating file names. This construction
# is useful when files may be moved to different directories.
set(module_list Matrix
  Vector
  BasisReader
  BasisWriter
  SVDBasisGenerator
  SVDSampler
  SVD
  IncrementalSVD
  IncrementalSVDFastUpdate
  IncrementalSVDStandard
  IncrementalSVDBasisGenerator
  IncrementalSVDSampler
  StaticSVD
  StaticSVDBasisGenerator
  StaticSVDSampler
  Database
  HDFDatabase
  Utilities
  ParallelBuffer
  DEIM
  QDEIM
  GNAT)
set(source_files)
foreach(module IN LISTS module_list)
  list(APPEND source_files ${module}.C)
endforeach(module) # IN LISTS module_list

add_library(ROM ${source_files})

# PUBLIC dependencies are transitive; these dependencies are used in
# the API headers *and* in their implementations
#
# INTERFACE dependencies are used in the API headers, but not in the
# API implementation (e.g., API forwards objects to another library
# without dereferencing those objects -- only pointers are used)
#
# PRIVATE dependencies are used in the API implementation, but not in
# the headers
target_link_libraries(ROM
  PUBLIC MPI HDF5 LAPACK
  PRIVATE ZLIB::ZLIB)

target_include_directories(ROM PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${HDF5_C_INCLUDE_DIRS}
  ${MPI_C_INCLUDE_DIRS}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Use the C++11 standard as an entire feature instead of
# enumerating individual compiler features for simplicity
target_compile_features(ROM PRIVATE cxx_std_11)

set(regression_test_names
  smoke_test
  uneven_dist
  weak_scaling
  random_test)

foreach(name IN LISTS regression_test_names)
  add_executable(${name} ${name}.C)
  target_link_libraries(${name} PRIVATE ROM MPI)
  target_include_directories(${name}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${MPI_C_INCLUDE_DIRS})
  target_compile_features(${name} PRIVATE cxx_std_11)
endforeach(name) # IN LISTS regression_test_names

if(GTEST_FOUND)
  set(unit_test_stems
    Vector
    Matrix
    SVD
    IncrementalSVD)
  foreach(stem IN LISTS unit_test_stems)
    add_executable(test_${stem} tests/test_${stem}.C)
    target_link_libraries(test_${stem} PRIVATE ROM MPI GTest::GTest)
    target_include_directories(test_${stem}
      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
      ${MPI_C_INCLUDE_DIRS})
    target_compile_features(test_${stem} PRIVATE cxx_std_11)
    target_compile_definitions(test_${stem} PRIVATE CAROM_HAS_GTEST)
  endforeach(stem) # IN LISTS unit_test_stems
endif(GTEST_FOUND)

# NOTE(goxberry@gmail.com, oxberry1@llnl.gov): This code snippet
# builds the Doxygen documentation, but outputs said documentation to
# ${CMAKE_CURRENT_SOURCE_DIR}/docs/rom-dox)...
if(DOXYGEN_FOUND)
  set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
  set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

  configure_file(${doxyfile_in} ${doxyfile} @ONLY)

  add_custom_target(
    documentation ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/docs
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM
    )
  add_dependencies(documentation ROM)

  add_custom_target(
    doxygen_tagfile
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/rom.tag
    ${CMAKE_CURRENT_BINARY_DIR}/rom-dox/html/rom.tag)
  add_dependencies(doxygen_tagfile documentation)

endif(DOXYGEN_FOUND)
