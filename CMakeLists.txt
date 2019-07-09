# Require CMake version 3.10 or later to expose modern features of
# CMake for maintainability. With additional work, it is possible to
# rewrite the build system generator for CMake 3.8+, 3.1+, 3.0+, or
# 2.8+, with increasing amounts of work as the minimum required
# version is decreased.
cmake_minimum_required(VERSION 3.8)

project(librom
  VERSION 0.0.1
  DESCRIPTION "Model reduction library emphasizing large-scale parallelism"
  HOMEPAGE_URL "https://github.com/LLNL/libROM"
  LANGUAGES C CXX)

find_package(HDF5 1.8.0
  REQUIRED)

find_package(LAPACK 3.4.0
  REQUIRED)

# In FindMPI.cmake, "CXX" refers to "the MPI C API [being] usable from C++"
find_package(MPI 1.2 REQUIRED)

find_package(ZLIB 1.2.3 REQUIRED)

find_package(Doxygen 1.8.10)

# find_package(Perl 5.16.3 REQUIRED)

find_package(GTest 1.6.0)

add_library(libROM Matrix.C
  Vector.C
  BasisReader.C
  BasisWriter.C
  SVDBasisGenerator.C
  SVDSampler.C
  SVD.C
  IncrementalSVD.C
  IncrementalSVDFastUpdate.C
  IncrementalSVDStandard.C
  IncrementalSVDBasisGenerator.C
  IncrementalSVDSampler.C
  StaticSVD.C
  StaticSVDBasisGenerator.C
  StaticSVDSampler.C
  Database.C
  HDFDatabase.C
  Utilities.C
  ParallelBuffer.C
  DEIM.C
  QDEIM.C
  GNAT.C)

# PUBLIC dependencies are transitive; these dependencies are used in
# the API headers *and* in their implementations
#
# INTERFACE dependencies are used in the API headers, but not in the
# API implementation (e.g., API forwards objects to another library
# without dereferencing those objects -- only pointers are used)
#
# PRIVATE dependencies are used in the API implementation, but not in
# the headers
target_link_libraries(libROM
  PUBLIC MPI HDF5 LAPACK
  PRIVATE ZLIB::ZLIB)

target_include_directories(libROM PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${HDF5_C_INCLUDE_DIRS}
  ${MPI_C_INCLUDE_DIRS}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Use the C++11 standard as an entire feature instead of
# enumerating individual compiler features for simplicity
target_compile_features(libROM PRIVATE cxx_std_11)

add_executable(smoke_test smoke_test.C)
target_link_libraries(smoke_test PRIVATE libROM MPI)
target_include_directories(smoke_test
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
  ${MPI_C_INCLUDE_DIRS})
target_compile_features(smoke_test PRIVATE cxx_std_11)

add_executable(uneven_dist uneven_dist.C)
target_link_libraries(uneven_dist PRIVATE libROM MPI)
target_include_directories(uneven_dist
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
  ${MPI_C_INCLUDE_DIRS})
target_compile_features(uneven_dist PRIVATE cxx_std_11)

add_executable(weak_scaling weak_scaling.C)
target_link_libraries(weak_scaling PRIVATE libROM MPI)
target_include_directories(weak_scaling
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
  ${MPI_C_INCLUDE_DIRS})
target_compile_features(weak_scaling PRIVATE cxx_std_11)

add_executable(random_test random_test.C)
target_link_libraries(random_test PRIVATE libROM MPI)
target_include_directories(random_test
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
  ${MPI_C_INCLUDE_DIRS})
target_compile_features(random_test PRIVATE cxx_std_11)

if(GTEST_FOUND)
  add_executable(test_Vector tests/test_Vector.C)
  target_link_libraries(test_Vector PRIVATE libROM MPI GTest::GTest)
  target_include_directories(test_Vector
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${MPI_C_INCLUDE_DIRS})
  target_compile_features(test_Vector PRIVATE cxx_std_11)
  target_compile_definitions(test_Vector PRIVATE CAROM_HAS_GTEST)

  add_executable(test_Matrix tests/test_Matrix.C)
  target_link_libraries(test_Matrix PRIVATE libROM MPI GTest::GTest)
  target_include_directories(test_Matrix
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${MPI_C_INCLUDE_DIRS})
  target_compile_features(test_Matrix PRIVATE cxx_std_11)
  target_compile_definitions(test_Matrix PRIVATE CAROM_HAS_GTEST)

  add_executable(test_SVD tests/test_SVD.C)
  target_link_libraries(test_SVD PRIVATE libROM MPI GTest::GTest)
  target_include_directories(test_SVD
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${MPI_C_INCLUDE_DIRS})
  target_compile_features(test_SVD PRIVATE cxx_std_11)
  target_compile_definitions(test_SVD PRIVATE CAROM_HAS_GTEST)

  add_executable(test_IncrementalSVD tests/test_IncrementalSVD.C)
  target_link_libraries(test_IncrementalSVD PRIVATE libROM MPI GTest::GTest)
  target_include_directories(test_IncrementalSVD
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${MPI_C_INCLUDE_DIRS})
  target_compile_features(test_IncrementalSVD PRIVATE cxx_std_11)
  target_compile_definitions(test_IncrementalSVD PRIVATE CAROM_HAS_GTEST)
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
  add_dependencies(documentation libROM)

  add_custom_target(
    doxygen_tagfile
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/rom.tag
    ${CMAKE_CURRENT_BINARY_DIR}/rom-dox/html/rom.tag)
  add_dependencies(doxygen_tagfile documentation)

endif(DOXYGEN_FOUND)
