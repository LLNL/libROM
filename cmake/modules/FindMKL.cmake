#.rst:
# FindMKL
# ---------
#
# Find Intel Math Kernel Library include dirs and libraries
#
# Use this module by invoking find_package with the form::
#
#   find_package(MKL
#     [version] [EXACT]      # Minimum or EXACT version e.g. 11.3.2,
#     [REQUIRED]             # Fail with error if MKL is not found
#     [COMPONENTS <libs>...] # Cluster libraries e.g. "PARDISO", "CDFT",
#     )                      # "ScaLAPACK" or "BLACS" (case sensitive).
#                            # Dependencies (PARDISO, CDFT and ScaLAPACK require
#                            # BLACS) are handled automatically
#
# This module finds headers and requested cluster libraries.  Results are
# reported in variables::
#
#   MKL_FOUND            - True if headers and requested libraries were found
#   MKL_DEFINITIONS      - MKL compiler definitions
#   MKL_INCLUDE_DIRS     - MKL include directories
#   MKL_LIBRARIES        - MKL libraries to be linked
#   MKL_<C>_FOUND        - True if cluster library <C> was found (<C> is
#                          upper-case)
#   MKL_<C>_LIBRARY      - Libraries to link for cluster library <C> (may
#                          include target_link_libraries debug/optimized
#                          keywords)
#   MKL_VERSION          - INTEL_MKL_VERSION value from mkl.h
#   MKL_MAJOR_VERSION    - MKL major version number (X in X.y.z)
#   MKL_MINOR_VERSION    - MKL minor version number (Y in x.Y.z)
#   MKL_PATCH_VERSION    - MKL patch version number (Z in x.y.Z)
#
# This module reads hints about search locations from variables::
#
#   MKL_ROOT             - Preferred installation prefix
#    (or MKLROOT)
#   MKL_INCLUDEDIR       - Preferred include directory e.g. <prefix>/include
#   MKL_LIBRARYDIR       - Preferred library directory e.g. <prefix>/lib
#   MKL_ADDITIONAL_VERSIONS -
#
# and saves search results persistently in CMake cache entries::
#
#   MKL_INCLUDE_DIR         - Directory containing MKL headers
#   MKL_LIBRARY_DIR         - Directory containing MKL libraries
#   MKL_<C>_LIBRARY         - Component <C> library

# This module first searches for mkl.h using the above hint variables (excluding
# MKL_LIBRARYDIR) and saves the result in MKL_INCLUDE_DIR.  Then it searches for
# appropriate interface libraries for the current architecture using dynamic or
# static linking.  Finally it searches for requested cluster libraries using
# the above hints (excluding MKL_INCLUDEDIR and MKL_ADDITIONAL_VERSIONS), "lib"
# directories near MKL_INCLUDE_DIR, and the library name configuration settings
# below.  It saves the library directories in MKL_LIBRARY_DIR and individual
# library locations in MKL_<C>_LIBRARY.
# When one changes settings used by previous searches in the same build
# tree (excluding environment variables) this module discards previous
# search results affected by the changes and searches again.
#
# MKL libraries come in many variants encoded in their file name.
# Users or projects may tell this module which variant to find by
# setting variables::
#
#   MKL_USE_STATIC_LIBS      - Set to ON to force the use of static libraries.
#                              Default is OFF.
#   MKL_INTERFACE            - Set to "lp64" or "ilp64" (case-insensitive) to
#                              select 32-bit or 64-bit integer interface.
#                              Default is blank on 32-bit platforms (where no
#                              interface library is needed) and to guess based
#                              on sizeof(int) on 64-bit platforms.
#   MKL_THREADING            - Set to "Sequential", "OpenMP" or "TBB"  to select
#                              threading library.  Default is "OpenMP".
#   MKL_OPENMP               - Set to "Intel", "GNU" or "PGI" to use vendor
#                              OpenMP library.  This setting only has an effect
#                              when MPI_THREADING is"OpenMP".  Default is
#                              "Intel".
#   MKL_MESSAGE_PASSING      - Set to "Intel", "MPICH", "MPICH2", "OpenMPI" or
#                              "SGI" to use a specific MPI library.  This
#                              setting only has an effect when COMPONENTS is
#                              non-empty.  Default is "Intel".
#
# Some combinations of the variants above do not exist (e.g. Using GNU compilers
# with PGI OpenMP and SGI Message Passing).  In order to support future versions
# of MKL these combinations are not checked for in CMake and will result in a
# generic "library not found" error.
#
# Example to find MKL headers only::
#
#   find_package(MKL 11.3.2)
#   if(MKL_FOUND)
#     include_directories(${MKL_INCLUDE_DIRS})
#     add_executable(foo foo.cc)
#   endif()
#
#
# Example to find MKL headers and some *static* libraries::
#
#   set(MKL_USE_STATIC_LIBS        ON) # only find static libs
#   find_package(MKL 11.3.2 COMPONENTS CDFT ScaLAPACK ...)
#   if(MKL_FOUND)
#     include_directories(${MKL_INCLUDE_DIRS})
#     link_directories(${MKL_LIBRARY_DIRS})
#     add_executable(foo foo.cc)
#     target_link_libraries(foo ${MKL_LIBRARIES})
#   endif()
#

# Set defaults
if (NOT DEFINED MKL_USE_STATIC_LIBS)
  set(MKL_USE_STATIC_LIBS OFF)
endif()
if (NOT DEFINED MKL_THREADING)
  set(MKL_THREADING OpenMP)
endif()
if (NOT DEFINED MKL_OPENMP)
  set(MKL_OPENMP Intel)
endif()
if (NOT DEFINED MKL_MESSAGE_PASSING)
  set(MKL_MESSAGE_PASSING Intel)
endif()

# Detect changes in used variables.
# Compares the current variable value with the last one.
# In short form:
# v != v_LAST                      -> CHANGED = 1
# v is defined, v_LAST not         -> CHANGED = 1
# v is not defined, but v_LAST is  -> CHANGED = 1
# otherwise                        -> CHANGED = 0
# CHANGED is returned in variable named ${changed_var}
macro(_MKL_CHANGE_DETECT changed_var)
  set(${changed_var} 0)
  foreach(v ${ARGN})
    if(DEFINED _MKL_COMPONENTS_SEARCHED)
      if(${v})
        if(_${v}_LAST)
          string(COMPARE NOTEQUAL "${${v}}" "${_${v}_LAST}" _${v}_CHANGED)
        else()
          set(_${v}_CHANGED 1)
        endif()
      elseif(_${v}_LAST)
        set(_${v}_CHANGED 1)
      endif()
      if(_${v}_CHANGED)
        set(${changed_var} 1)
      endif()
    else()
      set(_${v}_CHANGED 0)
    endif()
  endforeach()
endmacro()

# Collect environment variable inputs as hints.  Do not consider changes.
foreach(v MKLROOT MKL_ROOT MKL_INCLUDEDIR MKL_LIBRARYDIR)
  set(_env $ENV{${v}})
  if(_env)
    file(TO_CMAKE_PATH "${_env}" _ENV_${v})
  else()
    set(_ENV_${v} "")
  endif()
endforeach()
if(NOT _ENV_MKL_ROOT AND _ENV_MKLROOT)
  set(_ENV_MKL_ROOT "${_ENV_MKLROOT}")
endif()

# Collect inputs and cached results.  Detect changes since the last run.
if(NOT MKL_ROOT AND MKLROOT)
  set(MKL_ROOT "${MKLROOT}")
endif()
set(_MKL_VARS_DIR
  MKL_ROOT
  MKL_NO_SYSTEM_PATHS
  )

# ------------------------------------------------------------------------
#  Search for MKL include DIR
# ------------------------------------------------------------------------

set(_MKL_VARS_INC MKL_INCLUDEDIR MKL_INCLUDE_DIR MKL_ADDITIONAL_VERSIONS)
_MKL_CHANGE_DETECT(_MKL_CHANGE_INCDIR ${_MKL_VARS_DIR} ${_MKL_VARS_INC})
# Clear MKL_INCLUDE_DIR if it did not change but other input affecting the
# location did.  We will find a new one based on the new inputs.
if(_MKL_CHANGE_INCDIR AND NOT _MKL_INCLUDE_DIR_CHANGED)
  unset(MKL_INCLUDE_DIR CACHE)
endif()

if(NOT MKL_INCLUDE_DIR)
  set(_MKL_INCLUDE_SEARCH_DIRS "")
  if(MKL_INCLUDEDIR)
    list(APPEND _MKL_INCLUDE_SEARCH_DIRS ${MKL_INCLUDEDIR})
  elseif(_ENV_MKL_INCLUDEDIR)
    list(APPEND _MKL_INCLUDE_SEARCH_DIRS ${_ENV_MKL_INCLUDEDIR})
  endif()

  if(MKL_ROOT)
    list(APPEND _MKL_INCLUDE_SEARCH_DIRS ${MKL_ROOT}/include ${MKL_ROOT})
  elseif(_ENV_MKL_ROOT)
    list(APPEND _MKL_INCLUDE_SEARCH_DIRS ${_ENV_MKL_ROOT}/include ${_ENV_MKL_ROOT})
  endif()

  if(MKL_NO_SYSTEM_PATHS)
    list(APPEND _MKL_INCLUDE_SEARCH_DIRS NO_CMAKE_SYSTEM_PATH)
  else()
    list(APPEND _MKL_INCLUDE_SEARCH_DIRS PATHS
      "C:/Program Files (x86)/IntelSWTools"
      /opt/intel
      )
  endif()

  # Create a list of directories to search.
  # MKL comes bundled with Intel Composer/Parallel Studio and has a different
  # version to that of the package it comes bundled with (e.g. Composer XE
  # 2016.1.111 bundles MKL 11.3.2)
  # TODO: Prepend new versions to this as they are released (they must remain
  # sorted in descending order).
  set(_MKL_PACKAGE_KNOWN_VERSIONS ${MKL_PACKAGE_ADDITIONAL_VERSIONS}
      "2016.2.181" "2016.1.111" "2016.0.109" "2015.2.164")
  foreach(v ${_MKL_PACKAGE_KNOWN_VERSIONS})
    if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
      list(APPEND _MKL_PATH_SUFFIXES "compilers_and_libraries_${v}/windows/mkl/include")
    elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      list(APPEND _MKL_PATH_SUFFIXES "compilers_and_libraries_${v}/mac/mkl/include")
    elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
      list(APPEND _MKL_PATH_SUFFIXES "compilers_and_libraries_${v}/linux/mkl/include")
    endif()
    list(APPEND _MKL_PATH_SUFFIXES "composer_xe_${v}/mkl/include")
  endforeach()

  find_path(MKL_INCLUDE_DIR
    NAMES         mkl.h
    HINTS         ${_MKL_INCLUDE_SEARCH_DIRS}
    PATH_SUFFIXES ${_MKL_PATH_SUFFIXES}
  )
endif()

# ------------------------------------------------------------------------------
#  Extract version information from mkl.h (or mkl_version.h in newer versions)
# ------------------------------------------------------------------------------

# Set MKL_FOUND based only on header location and version.
# It will be updated below for component libraries.
if (MKL_INCLUDE_DIR)
  # Extract __INTEL_MKL__, __INTEL_MKL_MINOR__ and __INTEL_MKL_UPDATE__ from
  # mkl.h or mkl_version.h
  if(EXISTS "${MKL_INCLUDE_DIR}/mkl_version.h")
    file(STRINGS "${MKL_INCLUDE_DIR}/mkl_version.h" _MKL_VERSION_CONTENTS REGEX "#define __INTEL_MKL(_MINOR|_UPDATE)?__ ")
  else()
    file(STRINGS "${MKL_INCLUDE_DIR}/mkl.h" _MKL_VERSION_CONTENTS REGEX "#define __INTEL_MKL(_MINOR|_UPDATE)?__ ")
  endif()

  if("${_MKL_VERSION_CONTENTS}" MATCHES "#define __INTEL_MKL__[ \t\r\n]+([0-9]+)")
    set(MKL_MAJOR_VERSION "${CMAKE_MATCH_1}")
  endif()
  if("${_MKL_VERSION_CONTENTS}" MATCHES "#define __INTEL_MKL_MINOR__[ \t\r\n]+([0-9]+)")
    set(MKL_MINOR_VERSION "${CMAKE_MATCH_1}")
  endif()
  if("${_MKL_VERSION_CONTENTS}" MATCHES "#define __INTEL_MKL_UPDATE__[ \t\r\n]+([0-9]+)")
    set(MKL_UPDATE_VERSION "${CMAKE_MATCH_1}")
  endif()
  unset(_MKL_VERSION_CONTENTS)

  math(EXPR MKL_VERSION "${MKL_MAJOR_VERSION} * 10000 + ${MKL_MINOR_VERSION} * 100 + ${MKL_UPDATE_VERSION}")
  set(_MKL_VERSION "${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}.${MKL_UPDATE_VERSION}")
  if (NOT MKL_FIND_QUIETLY)
    message(STATUS "Found MKL: ${MKL_INCLUDE_DIR} (found version \"${_MKL_VERSION}\")")
  endif()

  # Check version against any requested
  if (MKL_FIND_VERSION_EXACT AND NOT "${_MKL_VERSION}" VERSION_EQUAL "${MKL_FIND_VERSION}")
    set(MKL_FOUND 0)
  elseif (MKL_FIND_VERSION AND "${_MKL_VERSION}" VERSION_LESS "${MKL_FIND_VERSION}")
    set(MKL_FOUND 0)
  else()
    set(MKL_FOUND 1)
    set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
  endif()
  unset(_MKL_VERSION)
endif()

# ------------------------------------------------------------------------------
# Find libraries
# ------------------------------------------------------------------------------
if (MKL_FOUND)

  # Work out whether to search the ia32/ or intel64/ lib/ subdirectories
  set(_MKL_LIBRARY_SEARCH_DIRS
      "${MKL_INCLUDE_DIR}/../lib;${MKL_INCLUDE_DIR}/../../compiler/lib;${MKL_INCLUDE_DIR}/../../tbb/lib")
  try_run(_MKL_IS_64BIT
          _MKL_IS_64BIT_COMPILE_RESULT
          "${CMAKE_BINARY_DIR}"
          "${CMAKE_CURRENT_LIST_DIR}/arch/test_is_64bit.c")
  if (_MKL_IS_64BIT)
    list(APPEND _MKL_LIBRARY_SEARCH_DIRS "${MKL_INCLUDE_DIR}/../lib/intel64;${MKL_INCLUDE_DIR}/../../compiler/lib/intel64;${MKL_INCLUDE_DIR}/../../tbb/lib/intel64")
  else()
    list(APPEND _MKL_LIBRARY_SEARCH_DIRS "${MKL_INCLUDE_DIR}/../lib/ia32;${MKL_INCLUDE_DIR}/../../compiler/lib/ia32;${MKL_INCLUDE_DIR}/../../tbb/lib/intel64")
  endif()

  set(MKL_LIBRARIES "")
  set(MKL_DEFINITIONS "")

  # Handle static libraries
  set(_MKL_LIB_PREFIX "")
  set(_MKL_LIB_SUFFIX "")
  if(MKL_USE_STATIC_LIBS)
    set(_MKL_LIB_PREFIX "lib")
    set(_MKL_LIB_SUFFIX ".a")
    list(APPEND MKL_LIBRARIES "-Wl,--start-group")
  endif()

  # Find the core library
  find_library(MKL_CORE_LIB
               "${_MKL_LIB_PREFIX}mkl_core${_MKL_LIB_SUFFIX}"
               HINTS ${_MKL_LIBRARY_SEARCH_DIRS})
  if("${MKL_CORE_LIB}" STREQUAL "MKL_CORE_LIB-NOTFOUND")
    set(MKL_FOUND 0)
  else()
    list(APPEND MKL_LIBRARIES ${MKL_CORE_LIB})
  endif()

  # Work out interface layer
  if (NOT DEFINED MKL_INTERFACE)
    if (_MKL_IS_64BIT)
      include(CheckTypeSize)
      CHECK_TYPE_SIZE("int" _MKL_SIZEOF_INT)
      if (${_MKL_SIZEOF_INT} EQUAL 4)
        set(MKL_INTERFACE "lp64")
      elseif(${_MKL_SIZEOF_INT} EQUAL 8)
        set(MKL_INTERFACE "ilp64")
      endif()
      unset(_MKL_SIZEOF_INT)
    endif()
  endif()
  if (DEFINED MKL_INTERFACE)
    string(TOLOWER ${MKL_INTERFACE} _MKL_INTERFACE_LC)
    if (${CMAKE_COMPILER_IS_GNUG77})
      set(_MKL_INTERFACE_LIBRARY "mkl_gf_${_MKL_INTERFACE_LC}")
    else()
      set(_MKL_INTERFACE_LIBRARY "mkl_intel_${_MKL_INTERFACE_LC}")
    endif()
    find_library(MKL_INTERFACE_LIB
                 "${_MKL_LIB_PREFIX}${_MKL_INTERFACE_LIBRARY}${_MKL_LIB_SUFFIX}"
                 HINTS ${_MKL_LIBRARY_SEARCH_DIRS})
    if("${MKL_INTERFACE_LIB}" STREQUAL "MKL_INTERFACE_LIB-NOTFOUND")
      set(MKL_FOUND 0)
    else()
      if("${_MKL_INTERFACE_LC}" STREQUAL "ilp64")
        list(APPEND MKL_DEFINITIONS "-DMKL_ILP64")
      endif()
      list(APPEND MKL_LIBRARIES ${MKL_INTERFACE_LIB})
    endif()
    unset(_MKL_INTERFACE_LIBRARY)
    unset(_MKL_IS_64BIT)
  endif()

  # Work out threading layer
  cmake_policy(PUSH)
  cmake_policy(SET CMP0054 NEW)
  string(TOLOWER ${MKL_THREADING} _MKL_THREADING_LC)
  if(_MKL_THREADING_LC STREQUAL "sequential")
    set(_MKL_THREADING_LIBS "mkl_sequential")
  elseif(_MKL_THREADING_LC STREQUAL "tbb")
    set(_MKL_THREADING_LIBS "mkl_tbb_thread;tbb")
  elseif(_MKL_THREADING_LC STREQUAL "openmp")
    string(TOLOWER ${MKL_OPENMP} _MKL_OPENMP_LC)
    if (_MKL_OPENMP_LC STREQUAL "intel")
      set(_MKL_THREADING_LIBS "mkl_intel_thread;iomp5")
    elseif(_MKL_OPENMP_LC STREQUAL "gnu")
      set(_MKL_THREADING_LIBS "mkl_gnu_thread;gomp")
    elseif(_MKL_OPENMP_LC STREQUAL "pgi")
      set(_MKL_THREADING_LIBS "mkl_pgi_thread;pgf90libs")
    endif()
  endif()
  cmake_policy(POP)
  foreach(lib ${_MKL_THREADING_LIBS})
    find_library(_MKL_THREADING_${lib}
                 "${_MKL_LIB_PREFIX}${lib}${_MKL_LIB_SUFFIX}"
                 HINTS ${_MKL_LIBRARY_SEARCH_DIRS})
    if("${_MKL_THREADING_${lib}}" STREQUAL "_MKL_THREADING_${lib}-NOTFOUND")
      set(MKL_FOUND 0)
    else()
      list(APPEND MKL_LIBRARIES ${_MKL_THREADING_${lib}})
    endif()
  endforeach()
  unset(_MKL_THREADING_LC)
  unset(_MKL_THREADING_LIBS)

  # All components currently require BLACS so if any are listed find BLACS first
  list(LENGTH MKL_FIND_COMPONENTS _MKL_NUM_COMPONENTS)
  if(${_MKL_NUM_COMPONENTS} GREATER 0)
    STRING(TOLOWER ${MKL_MESSAGE_PASSING} _MKL_MPI_LC)
    if(${_MKL_MPI_LC} STREQUAL "intel")
      set(_MKL_BLACS_LIBNAME "mkl_blacs_intelmpi_${_MKL_INTERFACE_LC}")
    elseif(${_MKL_MPI_LC} STREQUAL "mpich")
      set(_MKL_BLACS_LIBNAME "mkl_blacs_mpich_${_MKL_INTERFACE_LC}")
    elseif(${_MKL_MPI_LC} STREQUAL "mpich2")
      set(_MKL_BLACS_LIBNAME "mkl_blacs_intelmpi_${_MKL_INTERFACE_LC}")
    elseif(${_MKL_MPI_LC} STREQUAL "openmpi")
      set(_MKL_BLACS_LIBNAME "mkl_blacs_openmpi_${_MKL_INTERFACE_LC}")
    elseif(${_MKL_MPI_LC} STREQUAL "sgi")
      set(_MKL_BLACS_LIBNAME "mkl_blacs_sgimpt_${_MKL_INTERFACE_LC}")
    endif()
    unset(_MKL_MPI_LC)

    find_library(MKL_BLACS_LIBRARY
                 "${_MKL_LIB_PREFIX}${_MKL_BLACS_LIBNAME}${_MKL_LIB_SUFFIX}"
                 HINTS ${_MKL_LIBRARY_SEARCH_DIRS})
    if("${MKL_BLACS_LIBRARY}" STREQUAL "MKL_BLACS_LIBRARY-NOTFOUND")
      set(MKL_FOUND 0)
    else()
      list(APPEND MKL_LIBRARIES ${MKL_BLACS_LIBRARY})
      set(MKL_BLACS_FOUND TRUE)
    endif()

    # There is absolutely no pattern to the component libraries in MKL.
    # Cluster PARDISO has no library
    # CDFT is only available as libmkl_cdft_core.so
    # ScaLAPACK is available as libmkl_scalapack_lp64.so or
    # libmkl_scalapack_ilp64.so depending on MKL_INTERFACE and only on intel64
    # and not Mac OS

    # Cluster PARDISO
    list(FIND MKL_FIND_COMPONENTS PARDISO _MKL_INDEXOF_PARDISO)
    if (${_MKL_INDEXOF_PARDISO} GREATER -1)
      set(MKL_PARDISO_LIB "") # :trollface:
      set(MKL_PARDISO_FOUND TRUE)
    endif()
    unset(_MKL_INDEXOF_PARDISO)

    # CDFT
    list(FIND MKL_FIND_COMPONENTS CDFT _MKL_INDEXOF_CDFT)
    if (${_MKL_INDEXOF_CDFT} GREATER -1)
      find_library(MKL_CDFT_LIBRARY
                   "${_MKL_LIB_PREFIX}mkl_cdft_core${_MKL_LIB_SUFFIX}"
                   HINTS ${_MKL_LIBRARY_SEARCH_DIRS})
      if("${MKL_CDFT_LIBRARY}" STREQUAL "MKL_CDFT_LIBRARY-NOTFOUND")
        set(MKL_FOUND 0)
      else()
        list(APPEND MKL_LIBRARIES ${MKL_CDFT_LIBRARY})
        set(MKL_CDFT_FOUND TRUE)
      endif()
    endif()
    unset(_MKL_INDEXOF_CDFT)

    # ScaLAPACK
    list(FIND MKL_FIND_COMPONENTS ScaLAPACK _MKL_INDEXOF_SCALAPACK)
    if (${_MKL_INDEXOF_SCALAPACK} GREATER -1)
      find_library(MKL_ScaLAPACK_LIBRARY
                   "${_MKL_LIB_PREFIX}mkl_scalapack_${_MKL_INTERFACE_LC}${_MKL_LIB_SUFFIX}"
                   HINTS ${_MKL_LIBRARY_SEARCH_DIRS})
      if("${MKL_ScaLAPACK_LIBRARY}" STREQUAL "MKL_ScaLAPACK_LIBRARY-NOTFOUND")
        set(MKL_FOUND 0)
      else()
        list(APPEND MKL_LIBRARIES ${MKL_ScaLAPACK_LIBRARY})
        set(MKL_ScaLAPACK_FOUND TRUE)
      endif()
    endif()
    unset(_MKL_INDEXOF_SCALAPACK)
  endif()
  unset(_MKL_INTERFACE_LC)

  unset(_MKL_LIB_PREFIX)
  unset(_MKL_LIB_SUFFIX)

  if(MKL_USE_STATIC_LIBS)
    list(APPEND MKL_LIBRARIES "-Wl,--end-group")
  endif()

  # Construct MKL_LIBRARIES and MKL_DEFINITIONS
  list(REMOVE_DUPLICATES MKL_DEFINITIONS)
  list(REMOVE_DUPLICATES MKL_LIBRARIES)
endif()

# Add variables to cache
if(MKL_FOUND)
  set(MKL_DEFINITIONS  ${MKL_DEFINITIONS} CACHE STRING "MKL compiler definitions")
  set(MKL_INCLUDE_DIRS  ${MKL_INCLUDE_DIRS} CACHE STRING "MKL include directories")
  set(MKL_LIBRARIES   ${MKL_LIBRARIES} CACHE STRING "MKL libraries to be linked")
  if(MKL_PARDISO_FOUND)
    set(MKL_PARDISO_LIBRARY ${MKL_PARDISO_LIBRARY} CACHE STRING "PARDISO libraries")
  endif()
  if(MKL_CDFT_FOUND)
    set(MKL_CDFT_LIBRARY ${MKL_CDFT_LIBRARY} CACHE STRING "CDFT libraries")
  endif()
  if(MKL_ScaLAPACK_LIBRARY_FOUND)
    set(MKL_ScaLAPACK_LIBRARY ${MKL_ScaLAPACK_LIBRARY} CACHE STRING "ScaLAPACK libraries")
  endif()
  if(MKL_BLACS_LIBRARY_FOUND)
    set(MKL_BLACS_LIBRARY ${MKL_BLACS_LIBRARY} CACHE STRING "BLACS libraries")
  endif()
endif()
