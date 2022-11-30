# Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

# Function that converts a version string of the form 'major[.minor[.patch]]' to
# the integer ((major * 100) + minor) * 100 + patch.
function(librom_version_to_int VersionString VersionIntVar)
  if ("${VersionString}" MATCHES "^([0-9]+)(.*)$")
    set(Major "${CMAKE_MATCH_1}")
    set(MinorPatchString "${CMAKE_MATCH_2}")
  else()
    set(Major 0)
  endif()
  if ("${MinorPatchString}" MATCHES "^\\.([0-9]+)(.*)$")
    set(Minor "${CMAKE_MATCH_1}")
    set(PatchString "${CMAKE_MATCH_2}")
  else()
    set(Minor 0)
  endif()
  if ("${PatchString}" MATCHES "^\\.([0-9]+)(.*)$")
    set(Patch "${CMAKE_MATCH_1}")
  else()
    set(Patch 0)
  endif()
  math(EXPR VersionInt "(${Major}*100+${Minor})*100+${Patch}")
  set(${VersionIntVar} ${VersionInt} PARENT_SCOPE)
endfunction()
