#[[
This file is a part of ScalaWRAP.

Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
of the MIT license. See the top-level COPYRIGHT file for details.

SPDX-License-Identifier: MIT
]]

add_executable(all-tests test-main.cpp utils.cpp colmajor-global-colmajor.cpp
               colmajor-global-rowmajor.cpp rowmajor-global-colmajor.cpp
               rowmajor-global-rowmajor.cpp test-svd.cpp test-gemm.cpp)
target_include_directories(all-tests PRIVATE ${MPI_C_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(all-tests scalawrap)

