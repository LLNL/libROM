/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: A class implementing interface of SVD for the Randomized SVD
//              algorithm.

#include "RandomizedSVD.h"

#include "mpi.h"
#include "linalg/scalapack_wrapper.h"
#include "utils/mpi_utils.h"

#include <limits.h>
#include <algorithm>
#include <random>

#include <stdio.h>
#include <string.h>

namespace CAROM {

RandomizedSVD::RandomizedSVD(
    Options options) :
    StaticSVD(options),
    d_subspace_dim(options.randomized_subspace_dim) {
    srand(options.random_seed);
}

void
RandomizedSVD::computeSVD()
{
    d_samples->n = d_num_samples;
    delete_factorizer();

    const int num_rows = d_total_dim;
    const int num_cols = d_num_samples;
    if (d_subspace_dim < 1 || d_subspace_dim > std::min(num_rows, num_cols)) {
        d_subspace_dim = std::min(num_rows, num_cols);
    }

    // Get snapshot matrix in distributed format.
    // If there are less dimensions than samples, use the transpose instead.
    Matrix* snapshot_matrix;
    if (num_rows > num_cols) {
        snapshot_matrix = new Matrix(d_dims[static_cast<unsigned>(d_rank)],
                                     num_cols, true);
        for (int rank = 0; rank < d_num_procs; ++rank) {
            // gather_transposed_block does the same as gather_block, but transposes
            // it; here, it is used to go from column-major to row-major order.
            gather_transposed_block(&snapshot_matrix->item(0, 0), d_samples.get(),
                                    d_istarts[static_cast<unsigned>(rank)] + 1, 1,
                                    d_dims[static_cast<unsigned>(rank)], num_cols,
                                    rank);
        }
    }
    else {
        std::vector<int> snapshot_transpose_row_offset;
        const int num_transposed_rows = split_dimension(num_cols, MPI_COMM_WORLD);
        const int num_cols_check = get_global_offsets(num_transposed_rows,
                                   snapshot_transpose_row_offset,
                                   MPI_COMM_WORLD);
        CAROM_VERIFY(num_cols == num_cols_check);

        snapshot_matrix = new Matrix(num_transposed_rows,
                                     num_rows, true);

        for (int rank = 0; rank < d_num_procs; ++rank) {
            gather_block(&snapshot_matrix->item(0, 0), d_samples.get(),
                         1, snapshot_transpose_row_offset[rank] + 1,
                         num_rows, snapshot_transpose_row_offset[rank + 1] -
                         snapshot_transpose_row_offset[rank],
                         rank);
        }
    }

    int snapshot_matrix_distributed_rows = std::max(num_rows, num_cols);

    // Create a random matrix of smaller dimension to project the snapshot matrix
    // If debug mode is turned on, just set rand_mat as an identity matrix of smaller size
    // for reproducibility
    Matrix* rand_mat;
    if (d_debug_algorithm) {
        rand_mat = new Matrix(snapshot_matrix->numColumns(), d_subspace_dim, false);
        for (int i = 0; i < snapshot_matrix->numColumns(); i++) {
            for (int j = 0; j < d_subspace_dim; j++) {
                rand_mat->item(i, j) = (i == j);
            }
        }
    }
    else {
        rand_mat = new Matrix(snapshot_matrix->numColumns(), d_subspace_dim, false,
                              true);
    }

    // Project snapshot matrix onto random subspace
    Matrix* rand_proj = snapshot_matrix->mult(rand_mat);
    int rand_proj_rows = rand_proj->numRows();
    delete rand_mat;

    Matrix* Q = rand_proj->qr_factorize();

    // Project d_samples onto Q
    Matrix* svd_input_mat = Q->transposeMult(snapshot_matrix);
    int svd_input_mat_distributed_rows = svd_input_mat->numDistributedRows();

    SLPK_Matrix svd_input;
    initialize_matrix(&svd_input, svd_input_mat->numColumns(),
                      svd_input_mat->numDistributedRows(),
                      d_npcol, d_nprow, d_blocksize, d_blocksize);
    scatter_block(&svd_input, 1, 1,
                  svd_input_mat->getData(),
                  svd_input_mat->numColumns(), svd_input_mat->numDistributedRows(),0);
    delete svd_input_mat;
    delete snapshot_matrix;

    // This block does the actual ScaLAPACK call to do the factorization.
    svd_init(d_factorizer.get(), &svd_input);

    d_factorizer->dov = 1;

    factorize(d_factorizer.get());
    free_matrix_data(&svd_input);

    // Compute how many basis vectors we will actually return.
    int sigma_cutoff = 0, hard_cutoff = num_cols;
    if (d_singular_value_tol == 0) {
        sigma_cutoff = std::numeric_limits<int>::max();
    } else {
        for (int i = 0; i < num_cols; ++i) {
            if (d_factorizer->S[i] / d_factorizer->S[0] > d_singular_value_tol) {
                sigma_cutoff += 1;
            } else {
                break;
            }
        }
    }
    if (d_max_basis_dimension != -1 && d_max_basis_dimension < hard_cutoff) {
        hard_cutoff = d_max_basis_dimension;
    }
    int ncolumns = hard_cutoff < sigma_cutoff ? hard_cutoff : sigma_cutoff;
    CAROM_VERIFY(ncolumns >= 0);
    ncolumns = std::min(ncolumns, d_subspace_dim);

    // Allocate the appropriate matrices and gather their elements.
    d_S = new Vector(ncolumns, false);
    {
        CAROM_VERIFY(ncolumns >= 0);
        unsigned nc = static_cast<unsigned>(ncolumns);
        memset(&d_S->item(0), 0, nc*sizeof(double));
    }

    d_basis = new Matrix(svd_input_mat_distributed_rows, ncolumns, false);
    d_basis_right = new Matrix(d_subspace_dim, ncolumns, false);
    // Since the input to the SVD was transposed, U and V are switched.
    for (int rank = 0; rank < d_num_procs; ++rank) {
        // V is computed in the transposed order so no reordering necessary.
        gather_block(&d_basis->item(0, 0), d_factorizer->V, 1, 1,
                     ncolumns, d_subspace_dim, rank);
        // gather_transposed_block does the same as gather_block, but transposes
        // it; here, it is used to go from column-major to row-major order.
        gather_transposed_block(&d_basis_right->item(0, 0), d_factorizer->U,
                                1,
                                1, svd_input_mat_distributed_rows,
                                ncolumns, rank);
    }

    for (int i = 0; i < ncolumns; ++i)
        d_S->item(i) = d_factorizer->S[static_cast<unsigned>(i)];

    // Lift solution back to higher dimension
    Matrix* d_new_basis = Q->mult(d_basis);
    delete d_basis;
    d_basis = d_new_basis;

    if (num_rows <= num_cols) {
        Matrix* temp = d_basis;
        d_basis = d_basis_right;
        d_basis_right = temp;

        int local_rows = -1;
        if (d_basis->numRows() == d_total_dim)
            local_rows = d_dim;
        else {
            printf("WARNING: RandomizedSVD::d_basis has a subspace dimension %d, different from sample dimension %d\n",
                   d_basis->numRows(), num_rows);
            local_rows = split_dimension(d_basis->numRows());
        }
        d_basis->distribute(local_rows);
        d_basis_right->gather();
    }

    d_this_interval_basis_current = true;
    delete Q;
    release_context(&svd_input);

    if (d_debug_algorithm) {
        if (d_rank == 0) {
            printf("Distribution of sampler's A and U:\n");
        }
        print_debug_info(d_samples.get());
        MPI_Barrier(MPI_COMM_WORLD);

        if (d_rank == 0) {
            printf("Distribution of sampler's V:\n");
        }
        print_debug_info(d_factorizer->V);
        MPI_Barrier(MPI_COMM_WORLD);

        if (d_rank == 0) {
            printf("Computed singular values: ");
            for (int i = 0; i < ncolumns; ++i)
                printf("%8.4E  ", d_factorizer->S[i]);
            printf("\n");
        }
    }

}

}
