/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Implementation of the SnapshotDMD class.

#include "manifold_interp/PCHIPInterpolator.h"
#include "SnapshotDMD.h"
#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include "utils/CSVDatabase.h"
#include <algorithm>

namespace CAROM {

SnapshotDMD::~SnapshotDMD()
{
}

void SnapshotDMD::train(int k, const Matrix* W0, double linearity_tol)
{
    CAROM_VERIFY(d_snapshots.size() > 0);

    if(k >= d_snapshots.size())
    {
        interpolateToNSnapshots(k + 1);
    }

    const Matrix* f_snapshots = getSnapshotMatrix();
    CAROM_VERIFY(f_snapshots->numColumns() > 1);
    CAROM_VERIFY(k > 0 && k <= f_snapshots->numColumns() - 1);
    d_energy_fraction = -1.0;
    d_k = k;
    constructDMD(f_snapshots, d_rank, d_num_procs, W0, linearity_tol);

    delete f_snapshots;
}

void SnapshotDMD::train(double energy_fraction, const Matrix* W0,
                        double linearity_tol)
{
    DMD::train(energy_fraction,W0,linearity_tol);
}

void SnapshotDMD::interpolateToNSnapshots(int n)
{
    PCHIPInterpolator* interp = new PCHIPInterpolator();
    std::vector<Vector*> new_snapshots;
    std::vector<Vector*> new_times;

    interp->interpolate(d_sampled_times,d_snapshots,n,new_times,new_snapshots);
    d_snapshots = new_snapshots;
    d_sampled_times = new_times;
    d_dt = d_sampled_times[2]->getData()[0]-d_sampled_times[1]->getData()[0];
}

}
