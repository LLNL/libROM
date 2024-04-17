/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Implementation of the AdaptiveDMD algorithm.

#include "AdaptiveDMD.h"
#include "manifold_interp/VectorInterpolator.h"
#include "manifold_interp/SnapshotInterpolator.h"
#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include <algorithm>

namespace CAROM {

AdaptiveDMD::~AdaptiveDMD()
{
    for (auto interp_snapshot : d_interp_snapshots)
    {
        delete interp_snapshot;
    }
}

void AdaptiveDMD::train(double energy_fraction, const Matrix* W0,
                        double linearity_tol)
{
    const Matrix* f_snapshots = getInterpolatedSnapshots();
    CAROM_VERIFY(f_snapshots->numColumns() > 1);
    CAROM_VERIFY(energy_fraction > 0 && energy_fraction <= 1);
    d_energy_fraction = energy_fraction;
    constructDMD(f_snapshots, d_rank, d_num_procs, W0, linearity_tol);

    delete f_snapshots;
}

void AdaptiveDMD::train(int k, const Matrix* W0, double linearity_tol)
{
    const Matrix* f_snapshots = getInterpolatedSnapshots();
    CAROM_VERIFY(f_snapshots->numColumns() > 1);
    CAROM_VERIFY(k > 0 && k <= f_snapshots->numColumns() - 1);
    d_energy_fraction = -1.0;
    d_k = k;
    constructDMD(f_snapshots, d_rank, d_num_procs, W0, linearity_tol);

    delete f_snapshots;
}

void AdaptiveDMD::interpolateSnapshots()
{
    CAROM_VERIFY(d_interp_snapshots.size() == 0);
    CAROM_VERIFY(d_snapshots.size() == d_sampled_times.size());
    CAROM_VERIFY(d_sampled_times.size() > 1);

    SnapshotInterpolator* interp = new SnapshotInterpolator();
    std::vector<Vector*> new_times;

    d_interp_snapshots = interp->interpolate(d_sampled_times,d_snapshots,n,&new_times);
}

double AdaptiveDMD::getTrueDt() const
{
    return d_dt;
}

const Matrix* AdaptiveDMD::getInterpolatedSnapshots()
{
    if (d_interp_snapshots.size() == 0) interpolateSnapshots();
    return createSnapshotMatrix(d_interp_snapshots);
}

}
