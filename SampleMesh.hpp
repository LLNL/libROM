/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

#ifndef SAMPLEMESH_H
#define SAMPLEMESH_H

#include "mfem.hpp"
#include "mpi.h"

#include <set>

using namespace mfem;
using namespace std;

namespace CAROM {

    #define FULL_DOF_STENCIL

    void CreateSampleMesh(ParMesh& pmesh,
                          ParFiniteElementSpace& fespace1, ParFiniteElementSpace& fespace2,
                          const FiniteElementCollection& fecoll1, const FiniteElementCollection& fecoll2,
                          MPI_Comm& rom_com,
                          const vector<int>& sample_dofs,  /* block true DOFs for both spaces, for all processes */
                          const vector<int>& local_num_sample_dofs,
                          ParMesh*& sample_pmesh,
                          vector<int>& stencil_dofs, /* Local true DOF's on the original full mesh, restricted to the sample mesh stencil. */
                          vector<int>& all_stencil_dofs, /* stencil_dofs, gathered over all processes */
                          vector<int>& s2sp,
                          vector<int>& st2sp,
                          ParFiniteElementSpace*& spfespace1, ParFiniteElementSpace*& spfespace2);

  void GatherDistributedMatrixRows(const CAROM::Matrix& BR, const CAROM::Matrix& BW, const int rrdim, const int rwdim,
    #ifdef FULL_DOF_STENCIL
                                     const int NR, ParFiniteElementSpace& fespaceR, ParFiniteElementSpace& fespaceW,
    #endif
                                     const vector<int>& st2sp, const vector<int>& sprows,
				   const vector<int>& all_sprows, CAROM::Matrix& BRsp, CAROM::Matrix& BWsp);
}


#endif // SAMPLEMESH_H
