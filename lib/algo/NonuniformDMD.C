/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Implementation of the nonuniform DMD algorithm.

#include "NonuniformDMD.h"
#include "linalg/Matrix.h"

namespace CAROM {

NonuniformDMD::NonuniformDMD(int dim) : DMD(dim)
{}

std::pair<Matrix*, Matrix*>
NonuniformDMD::computePlusMinusSnapshotMatrices(const Matrix* snapshots)
{
    // TODO: Making two copies of the snapshot matrix has a lot of overhead.
    //       We need to figure out a way to do submatrix multiplication and to
    //       reimplement this algorithm using one snapshot matrix.
    Matrix* f_snapshots_minus = new Matrix(snapshots->numRows(),
                                           snapshots->numColumns() - 1, snapshots->distributed());
    Matrix* f_snapshots_plus = new Matrix(snapshots->numRows(),
                                          snapshots->numColumns() - 1, snapshots->distributed());

    // Break up snapshots into snapshots_minus and snapshots_plus
    // snapshots_minus = all columns of snapshots except last
    // snapshots_plus = all columns of snapshots except first
    for (int i = 0; i < snapshots->numRows(); i++)
    {
        for (int j = 0; j < snapshots->numColumns() - 1; j++)
        {
            f_snapshots_minus->item(i, j) = snapshots->item(i, j);
            f_snapshots_plus->item(i, j) = (snapshots->item(i, j + 1) - snapshots->item(i, j)) / (d_sampled_times[j + 1] ->item(0) - d_sampled_times[j]->item(0));
        }
    }

    return std::pair<Matrix*,Matrix*>(f_snapshots_minus, f_snapshots_plus);
}

void
NonuniformDMD::computePhi(struct DMDInternal dmd_internal_obj)
{
    d_phi_real = dmd_internal_obj.basis->mult(dmd_internal_obj.eigenpair->ev_real);
    d_phi_imaginary = dmd_internal_obj.basis->mult(dmd_internal_obj.eigenpair->ev_imaginary);
}

std::complex<double>
NonuniformDMD::computeEigExp(std::complex<double> eig, double t)
{
    return std::exp(t * eig);
}

}
