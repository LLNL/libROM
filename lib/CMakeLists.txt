###############################################################################
#
#  Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
#  and other libROM project developers. See the top-level COPYRIGHT
#  file for details.
#
#  SPDX-License-Identifier: (Apache-2.0 OR MIT)
#
###############################################################################

include(FortranCInterface)
FortranCInterface_HEADER(${CMAKE_CURRENT_SOURCE_DIR}/FCMangle.h
  MACRO_NAMESPACE "CAROM_FC_")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CAROM_config.h.in
  ${CMAKE_CURRENT_SOURCE_DIR}/CAROM_config.h @ONLY)

# While it is tempting to use file globbing here, file globbing is
# considered a "modern CMake anti-pattern" because "CMake is not a
# build system -- it is a build system _generator_". Instead, use
# some file name "stems" for generating file names. This construction
# is useful when files may be moved to different directories.
set(module_list
  linalg/BasisGenerator
  linalg/BasisReader
  linalg/BasisWriter
  linalg/Matrix
  linalg/Vector
  linalg/NNLS
  linalg/svd/IncrementalSVD
  linalg/svd/IncrementalSVDFastUpdate
  linalg/svd/IncrementalSVDStandard
  linalg/svd/RandomizedSVD
  linalg/svd/SVD
  linalg/svd/StaticSVD
  algo/DMD
  algo/AdaptiveDMD
  algo/NonuniformDMD
  algo/IncrementalDMD
  algo/DifferentialEvolution
  algo/greedy/GreedyCustomSampler
  algo/greedy/GreedyRandomSampler
  algo/greedy/GreedySampler
  algo/manifold_interp/Interpolator
  algo/manifold_interp/MatrixInterpolator
  algo/manifold_interp/VectorInterpolator
  hyperreduction/DEIM
  hyperreduction/GNAT
  hyperreduction/QDEIM
  hyperreduction/S_OPT
  hyperreduction/STSampling
  hyperreduction/Utilities
  utils/Database
  utils/HDFDatabase
  utils/CSVDatabase
  utils/Utilities
  utils/ParallelBuffer)
set(source_files)
foreach(module IN LISTS module_list)
  list(APPEND source_files ${module}.cpp)
endforeach(module) # IN LISTS module_list

list(APPEND source_files
  algo/ParametricDMD.h
  linalg/Options.h
  librom.h)

if (USE_MFEM)
    list(APPEND source_files
      mfem/PointwiseSnapshot.hpp
      mfem/PointwiseSnapshot.cpp
      mfem/Utilities.hpp
      mfem/Utilities.cpp
      mfem/SampleMesh.hpp
      mfem/SampleMesh.cpp)
endif()

list(APPEND source_files
  linalg/scalapack_c_wrapper.c
  linalg/scalapack_f_wrapper.f90)
if (BUILD_STATIC)
  add_library(ROM ${source_files})
else()
  add_library(ROM SHARED ${source_files})
endif()

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
