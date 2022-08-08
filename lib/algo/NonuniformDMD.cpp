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

NonuniformDMD::NonuniformDMD(int dim, Vector* state_offset, Vector* derivative_offset) : 
                             DMD(dim, state_offset, derivative_offset)
{}

NonuniformDMD::NonuniformDMD(std::string base_file_name) : DMD(base_file_name)
{}

NonuniformDMD::NonuniformDMD(std::vector<std::complex<double>> eigs,
                             Matrix* phi_real,
                             Matrix* phi_imaginary, int k,
                             bool in_offset, bool out_offset,
                             Vector* state_offset, Vector* derivative_offset,
                             double dt, double t_offset) : 
                             DMD(eigs, phi_real, phi_imaginary, k,
                                 in_offset, out_offset, state_offset, derivative_offset, dt, t_offset)
{}

std::pair<Matrix*, Matrix*>
NonuniformDMD::computeDMDSnapshotPair(const Matrix* snapshots)
{
    CAROM_VERIFY(snapshots->numColumns() > 1);

    // TODO: Making two copies of the snapshot matrix has a lot of overhead.
    //       We need to figure out a way to do submatrix multiplication and to
    //       reimplement this algorithm using one snapshot matrix.
    Matrix* f_snapshots_in = new Matrix(snapshots->numRows(),
                                        snapshots->numColumns() - 1, snapshots->distributed());
    Matrix* f_snapshots_out = new Matrix(snapshots->numRows(),
                                         snapshots->numColumns() - 1, snapshots->distributed());

    // Break up snapshots into snapshots_in and snapshots_out
    // snapshots_in = all columns of snapshots except last
    // snapshots_out = finite difference of all columns of snapshots
    for (int i = 0; i < snapshots->numRows(); i++)
    {
        for (int j = 0; j < snapshots->numColumns() - 1; j++)
        {
            f_snapshots_in->item(i, j) = snapshots->item(i, j);
            f_snapshots_out->item(i, j) = (snapshots->item(i, j + 1) - snapshots->item(i,j)) / 
                                          (d_sampled_times[j + 1]->item(0) - d_sampled_times[j]->item(0));
            if (d_in_offset) f_snapshots_in->item(i, j) -= d_state_offset->item(i);
            if (d_out_offset) f_snapshots_out->item(i, j) -= d_derivative_offset->item(i);
        }
    }

    return std::pair<Matrix*,Matrix*>(f_snapshots_in, f_snapshots_out);
}

void
NonuniformDMD::computePhi(struct DMDInternal dmd_internal_obj)
{
    d_phi_real = dmd_internal_obj.basis->mult(dmd_internal_obj.eigenpair->ev_real);
    d_phi_imaginary = dmd_internal_obj.basis->mult(
                          dmd_internal_obj.eigenpair->ev_imaginary);
}

std::complex<double>
NonuniformDMD::computeEigExp(std::complex<double> eig, double t)
{
    return std::exp(t * eig);
}

}
