/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Implementation of the nonuniform DMD algorithm.

#include "mfem.hpp"
#include "IncrementalDMD.h"
#include "linalg/Matrix.h"
#include "utils/Utilities.h"

using namespace mfem;

namespace CAROM {

IncrementalDMD::IncrementalDMD(int dim,
			       double dt,
			       Options svd_options,
			       std::string svd_base_file_name,
                               bool alt_output_basis,
                               Vector* state_offset) :
    DMD(dim, dt, alt_output_basis, state_offset)
{
    bg = new BasisGenerator(svd_options,
		    	    true,
		    	    svd_base_file_name);
}

IncrementalDMD::~IncrementalDMD()
{
    delete bg;
}

void
IncrementalDMD::load(std::string base_file_name)
{
    CAROM_ASSERT(!base_file_name.empty());

    DMD::load(base_file_name);
}

void
IncrementalDMD::save(std::string base_file_name)
{
    CAROM_ASSERT(!base_file_name.empty());
    CAROM_VERIFY(d_trained);

    DMD::save(base_file_name);
}

void
IncrementalDMD::train(int k, const Matrix* W0, double linearity_tol)
{
    const Matrix* f_snapshots = getSnapshotMatrix();
    CAROM_VERIFY(f_snapshots->numColumns() > 1);
    CAROM_VERIFY(k > 0 && k <= f_snapshots->numColumns() - 1);
    d_energy_fraction = -1.0;
    d_k = k;
    updateDMD(f_snapshots);
    delete f_snapshots;
}

void
IncrementalDMD::updateDMD(const Matrix* f_snapshots)
{
    std::pair<Matrix*, Matrix*> f_snapshot_pair = computeDMDSnapshotPair(f_snapshots);
    Matrix* f_snapshots_in = f_snapshot_pair.first;
    Matrix* f_snapshots_out = f_snapshot_pair.second;
    
    /* Incremental SVD
     * 
     * Does everything below all at once:
     * (1) Initialize SVD if not initialized
     * (2) Check linear dependence of new snapshot
     * (3) Rank-1 update
     *
    */
    double* u_in = f_snapshots_in->getColumn(f_snapshots_in->numColumns()-1)
	    			 ->getData();

    StopWatch svd_timer, rest_timer;

    svd_timer.Start();
    bg->takeSample(u_in, 0, d_dt, false); // what if norm(u_in) < eps at init?
    svd_timer.Stop();
    std::cout << "Time svd:" << svd_timer.RealTime() << std::endl;

    rest_timer.Start();
    bg->computeBasis();
    // Copy pointers from IncrementalSVD to DMD
    d_k = bg->getSingularValues()->dim(); // no. of singular values
    
    if (d_trained){ delete d_basis; } 
    Matrix* d_basis = new Matrix(*(bg->getSpatialBasis())); // left singular vectors
    std::vector<double> vec_values(bg->getSingularValues()->getData(),
		    		   bg->getSingularValues()->getData() + d_k);
    d_sv = vec_values; // singular values
    Matrix* d_basis_right = new Matrix(*(bg->getTemporalBasis())); // right singular vectors
    Matrix* d_S_inv = new Matrix(d_k, d_k, false);

    /*
     * Distributed vectors: skip for now
    
    for (int d_rank = 0; d_rank < d_num_procs; ++d_rank) {
        // V is computed in the transposed order so no reordering necessary.
        gather_block(&d_basis->item(0, 0), d_factorizer->V,
                     1, row_offset[static_cast<unsigned>(d_rank)]+1,
                     d_k, row_offset[static_cast<unsigned>(d_rank) + 1] -
                     row_offset[static_cast<unsigned>(d_rank)],
                     d_rank);

        // gather_transposed_block does the same as gather_block, but transposes
        // it; here, it is used to go from column-major to row-major order.
        gather_transposed_block(&d_basis_right->item(0, 0), d_factorizer->U, 1, 1,
                                f_snapshots_in->numColumns(), d_k, d_rank);
    }
    delete[] row_offset;

    */

    // Get inverse of singular values by multiplying by reciprocal.
    for (int i = 0; i < d_k; ++i)
    {
        d_S_inv->item(i, i) = 1 / d_sv[i];
    }

    Matrix* Q = NULL;

    /* Initial basis (W0): neglect for now
     *
    // If W0 is not null, we need to ensure it is in the basis of W.
    if (W0 != NULL)
    {
        CAROM_VERIFY(W0->numRows() == d_basis->numRows());
        CAROM_VERIFY(linearity_tol >= 0.0);
        std::vector<int> lin_independent_cols_W;

        // Copy W0 and orthogonalize.
        Matrix* d_basis_init = new Matrix(f_snapshots->numRows(), W0->numColumns(),
                                          true);
        for (int i = 0; i < d_basis_init->numRows(); i++)
        {
            for (int j = 0; j < W0->numColumns(); j++)
            {
                d_basis_init->item(i, j) = W0->item(i, j);
            }
        }
        d_basis_init->orthogonalize();

        Vector W_col(f_snapshots->numRows(), f_snapshots->distributed());
        Vector l(W0->numColumns(), true);
        Vector W0l(f_snapshots->numRows(), f_snapshots->distributed());
        // Find which columns of d_basis are linearly independent from W0
        for (int j = 0; j < d_basis->numColumns(); j++)
        {
            // l = W0' * u
            for (int i = 0; i < f_snapshots->numRows(); i++)
            {
                W_col.item(i) = d_basis->item(i, j);
            }
            d_basis_init->transposeMult(W_col, l);

            // W0l = W0 * l
            d_basis_init->mult(l, W0l);

            // Compute k = sqrt(u.u - 2.0*l.l + basisl.basisl) which is ||u -
            // basisl||_{2}.  This is the error in the projection of u into the
            // reduced order space and subsequent lifting back to the full
            // order space.
            double k = W_col.inner_product(W_col) - 2.0*l.inner_product(l) +
                       W0l.inner_product(W0l);
            if (k <= 0)
            {
                k = 0;
            }
            else
            {
         u.GetData       k = sqrt(k);
            }

            // Use k to see if the vector addressed by u is linearly dependent
            // on the left singular vectors.
            if (k >= linearity_tol)
            {
                lin_independent_cols_W.push_back(j);
            }
        }
        delete d_basis_init;

        // Add the linearly independent columns of W to W0. Call this new basis W_new.
        Matrix* d_basis_new = new Matrix(f_snapshots->numRows(),
                                         W0->numColumns() + lin_independent_cols_W.size(), true);
        for (int i = 0; i < d_basis_new->numRows(); i++)
        {
            for (int j = 0; j < W0->numColumns(); j++)
            {
                d_basis_new->item(i, j) = W0->item(i, j);
            }
            for (int j = 0; j < lin_independent_cols_W.size(); j++)
            {
                d_basis_new->item(i, j + W0->numColumns()) = d_basis->item(i,
                        lin_independent_cols_W[j]);
            }
        }

        // Orthogonalize W_new.
        d_basis_new->orthogonalize();

        // Calculate Q = W* x W_new;
        Q = d_basis->transposeMult(d_basis_new);

        delete d_basis;
        d_basis = d_basis_new;

        d_k = d_basis_new->numColumns();
        if (d_rank == 0) std::cout << "After adding W0, now using " << d_k <<
                                       " basis vectors." << std::endl;
    }

    */

    // Calculate A_tilde = U_transpose * f_snapshots_out * V * inv(S)
    Matrix* d_basis_mult_f_snapshots_out = d_basis->transposeMult(f_snapshots_out);
    Matrix* d_basis_mult_f_snapshots_out_mult_d_basis_right =
        d_basis_mult_f_snapshots_out->mult(d_basis_right);
    if (Q == NULL)
    {
        d_A_tilde = d_basis_mult_f_snapshots_out_mult_d_basis_right->mult(d_S_inv);
    }
    else
    {
        Matrix* d_basis_mult_f_snapshots_out_mult_d_basis_right_mult_d_S_inv =
            d_basis_mult_f_snapshots_out_mult_d_basis_right->mult(d_S_inv);
        d_A_tilde = d_basis_mult_f_snapshots_out_mult_d_basis_right_mult_d_S_inv->mult(
                        Q);
        delete Q;
        delete d_basis_mult_f_snapshots_out_mult_d_basis_right_mult_d_S_inv;
    }

    // Calculate the right eigenvalues/eigenvectors of A_tilde
    ComplexEigenPair eigenpair = NonSymmetricRightEigenSolve(d_A_tilde);
    d_eigs = eigenpair.eigs;

    struct DMDInternal dmd_internal = {f_snapshots_in, f_snapshots_out, d_basis, d_basis_right, d_S_inv, &eigenpair};
    computePhi(dmd_internal);

    Vector* init = new Vector(f_snapshots_in->numRows(), true);
    for (int i = 0; i < init->dim(); i++)
    {
        init->item(i) = f_snapshots_in->item(i, 0);
    }

    // Calculate pinv(d_phi) * initial_condition.
    projectInitialCondition(init);

    d_trained = true;

    delete d_basis_right;
    delete d_S_inv;
    delete d_basis_mult_f_snapshots_out;
    delete d_basis_mult_f_snapshots_out_mult_d_basis_right;
    delete f_snapshots_in;
    delete f_snapshots_out;
    delete eigenpair.ev_real;
    delete eigenpair.ev_imaginary;
    delete init;

    rest_timer.Stop();
    std::cout << "Time rest:" << rest_timer.RealTime() << std::endl;

}

}
