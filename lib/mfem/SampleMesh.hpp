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

#include <set>
#include <string>

#include "mfem.hpp"
#include "mpi.h"

#include "linalg/Matrix.h"

using namespace mfem;
using namespace std;

namespace CAROM {

#define FULL_DOF_STENCIL

class SampleMeshManager
{
public:
    SampleMeshManager(vector<ParFiniteElementSpace*> & fespace_, string visFileName);

    int RegisterSampledVariable(const int space, vector<int> const& sample_dofs_v, vector<int> const& num_sample_dofs_per_proc);
    void GatherDistributedMatrixRows(const int var, CAROM::Matrix const& B, const int rdim, CAROM::Matrix& Bsp);

    /**
     * @brief Destructor.
     */
    ~SampleMeshManager();

private:
    void ConstructSampleMesh();
    void CreateSampleMesh(vector<int>& stencil_dofs, vector<int>& all_stencil_dofs);

    void SetSampleMaps();
    void FinishSampleMaps();

    bool finalized;

    const int nspaces;
    int nvar;

    int myid, nprocs;
    MPI_Comm root_comm;

    vector<ParFiniteElementSpace*> fespace;
    vector<ParFiniteElementSpace*> spfespace;

    ParMesh *pmesh, *sample_pmesh;

    vector<vector<set<int>>> sample_dofs_proc;
    //vector<set<int>> sample_dofs_set;

    // TMP: sample_dofs is sample_dofs_merged in laghos_rom.cpp
    vector<int> sample_dofs;  // block true DOFs for all spaces, for all processes.

    vector<int> num_sample_dofs_per_proc_merged;

    vector<int> varSpace;
    vector<vector<int>> sample_dofs_var;
    vector<vector<int>> num_sample_dofs_per_proc_var;
    vector<vector<int>> s2sp_var;
    vector<int> s2sp, st2sp, sprows, all_sprows;

    vector<int> spaceTOS, spaceOS, spaceOSSP;
    vector<vector<int>> spaceOSall;

    vector<vector<int>> s2sp_space;

    // TODO: more descriptive name?
    string filename;  // For visualization output
};


}


#endif // SAMPLEMESH_H
