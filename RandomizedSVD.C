/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
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
#include "scalapack_wrapper.h"

#include <limits.h>
#include <algorithm>
#include <random>

#include <stdio.h>
#include <string.h>

/* Use automatically detected Fortran name-mangling scheme */
#define dormlq CAROM_FC_GLOBAL(dormlq, DORMLQ)

extern "C" {
void dormlq(char*, char*, int*, int*, int*,
    double*, int*, double*, double*, int*, double*, int*, int*);
}

namespace CAROM {

RandomizedSVD::RandomizedSVD(
   Options options) :
   StaticSVD(options),
   d_subspace_dim(options.randomized_subspace_dim) {
     srand(1);
   }

void
RandomizedSVD::computeSVD()
{
     d_samples->n = d_num_samples;
     delete_factorizer();

     int num_rows = d_total_dim;
     int num_cols = d_num_samples;
     if (d_subspace_dim < 1 || d_subspace_dim > std::min(num_rows, num_cols)) {
       d_subspace_dim = std::min(num_rows, num_cols);
     }

     // Get snapshot matrix in distributed format.
     // If there are less dimensions than samples, use the transpose instead.
     Matrix* snapshot_matrix;
     Matrix* snapshot_matrix_transposed;
     if (num_rows > num_cols) {
       snapshot_matrix = new Matrix(d_dims[static_cast<unsigned>(d_rank)],
           num_cols, true);
       snapshot_matrix_transposed = new Matrix(d_dims[static_cast<unsigned>(d_rank)],
          num_cols, true);
       for (int rank = 0; rank < d_num_procs; ++rank) {
          // gather_transposed_block does the same as gather_block, but transposes
          // it; here, it is used to go from column-major to row-major order.
          gather_transposed_block(&snapshot_matrix->item(0, 0), d_samples.get(),
                                 d_istarts[static_cast<unsigned>(rank)] + 1, 1,
                                 d_dims[static_cast<unsigned>(rank)], num_cols,
                                 rank);
          gather_block(&snapshot_matrix_transposed->item(0, 0), d_samples.get(),
                                 d_istarts[static_cast<unsigned>(rank)] + 1, 1,
                                 d_dims[static_cast<unsigned>(rank)], num_cols,
                                 rank);
       }
      }
      else {
        int num_transposed_rows = num_cols / d_num_procs;
        if (num_cols % d_num_procs > d_rank) {
          num_transposed_rows++;
        }
       snapshot_matrix = new Matrix(num_transposed_rows,
          num_rows, true);
       snapshot_matrix_transposed = new Matrix(num_transposed_rows,
          num_rows, true);
       for (int rank = 0; rank < d_num_procs; ++rank) {
         gather_block(&snapshot_matrix->item(0, 0), d_samples.get(),
                                d_istarts[static_cast<unsigned>(rank)] + 1, 1,
                                d_dims[static_cast<unsigned>(rank)], num_cols,
                                rank);
         gather_transposed_block(&snapshot_matrix_transposed->item(0, 0), d_samples.get(),
                                d_istarts[static_cast<unsigned>(rank)] + 1, 1,
                                d_dims[static_cast<unsigned>(rank)], num_cols,
                                rank);
       }
      }

     // Create a random matrix of smaller dimension to project the snapshot matrix
     Matrix* rand_mat = new Matrix(snapshot_matrix->numColumns(), d_subspace_dim, false, true);

     // Project snapshot matrix onto random subspace
     Matrix* rand_proj = snapshot_matrix->mult(rand_mat);

     // Get QR factorization of random projection
     int *row_offset = new int[d_num_procs + 1];
     row_offset[d_num_procs] = rand_proj->numDistributedRows();

     int my_rank;
     MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

     row_offset[my_rank] = rand_proj->numRows();

     CAROM_VERIFY(MPI_Allgather(MPI_IN_PLACE,
              1,
              MPI_INT,
              row_offset,
              1,
              MPI_INT,
              MPI_COMM_WORLD) == MPI_SUCCESS);
      for (int i = d_num_procs - 1; i >= 0; i--) {
        row_offset[i] = row_offset[i + 1] - row_offset[i];
      }

     CAROM_VERIFY(row_offset[0] == 0);

     SLPK_Matrix slpk_rand_proj;

     d_nprow = d_num_procs;
     d_npcol = 1;
     int d_blocksize = row_offset[d_num_procs] / d_num_procs;
     if (row_offset[d_num_procs] % d_num_procs != 0) d_blocksize += 1;
       initialize_matrix(&slpk_rand_proj, rand_proj->numDistributedRows(), rand_proj->numColumns(),
         d_nprow, d_npcol, d_blocksize, d_blocksize);
     for (int rank = 0; rank < d_num_procs; ++rank) {
        scatter_block(&slpk_rand_proj, row_offset[rank] + 1,
                  1, rand_proj->getData(), row_offset[rank + 1] - row_offset[rank],
                  rand_proj->numColumns(), rank);
     }

     QRManager QRmgr;
     qr_init(&QRmgr, &slpk_rand_proj);
     qrfactorize(&QRmgr);

     // Manipulate lq_A to get elementary household reflectors.
     Matrix* lq_A = new Matrix(QRmgr.A->mdata, d_subspace_dim,
        rand_proj->numRows(), false, false);

     for (int i = 0; i < lq_A->numRows(); i++) {
         for (int j = 0; j < i; j++) {
           lq_A->item(i, j) = 0;
         }
         lq_A->item(i, i) = 1;
     }

     SLPK_Matrix svd_input;
     make_similar_matrix(&svd_input, snapshot_matrix_transposed->numDistributedRows(),
            snapshot_matrix_transposed->numColumns(), QRmgr.A->ctxt, d_blocksize, d_blocksize);

     for (int rank = 0; rank < d_num_procs; ++rank) {
        scatter_block(&svd_input, 1, row_offset[rank] + 1,
                  snapshot_matrix_transposed->getData(), row_offset[rank + 1] - row_offset[rank],
                  snapshot_matrix_transposed->numColumns(), rank);

     }

    // This computes the action of Q on the input to the SVD.
    qaction(&QRmgr, &svd_input, 'L', 'N');

    // This block does the actual ScaLAPACK call to do the factorization.
    svd_init(d_factorizer.get(), &svd_input);

    d_factorizer->dov = 1;
    factorize(d_factorizer.get());

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
    d_basis = new Matrix(rand_proj->numRows(), ncolumns, true);
    d_S = new Matrix(ncolumns, ncolumns, false);
    {
       CAROM_VERIFY(ncolumns >= 0);
       unsigned nc = static_cast<unsigned>(ncolumns);
       memset(&d_S->item(0, 0), 0, nc*nc*sizeof(double));
    }
    d_basis_right = new Matrix(ncolumns, d_subspace_dim, false);

    SLPK_Matrix d_new_basis;
    make_similar_matrix(&d_new_basis, snapshot_matrix->numDistributedRows(),
           ncolumns, QRmgr.A->ctxt, d_blocksize, d_blocksize);
    for (int rank = 0; rank < d_num_procs; ++rank) {
       copy_matrix(&d_new_basis, row_offset[rank] + 1, 1,
                 d_factorizer->U, row_offset[rank] + 1, 1,
                 rand_proj->numRows(), ncolumns);
    }

    // This computes the action of Q on d_basis.
    qaction(&QRmgr, &d_new_basis, 'L', 'T');

    // Since the input to the SVD was transposed, U and V are switched.
    for (int rank = 0; rank < d_num_procs; ++rank) {
       // gather_transposed_block does the same as gather_block, but transposes
       // it; here, it is used to go from column-major to row-major order.
       gather_transposed_block(&d_basis->item(0, 0), &d_new_basis,
                               row_offset[rank] + 1,
                               1, rand_proj->numRows(),
                               ncolumns, rank);
       // V is computed in the transposed order so no reordering necessary.
       gather_block(&d_basis_right->item(0, 0), d_factorizer->V, 1, 1,
                    ncolumns, d_subspace_dim, rank);
    }
    for (int i = 0; i < ncolumns; ++i)
       d_S->item(i, i) = d_factorizer->S[static_cast<unsigned>(i)];
    if (num_rows <= num_cols) {
      Matrix* temp = d_basis;
      d_basis = d_basis_right;
      d_basis_right = temp;
    }
    d_this_interval_basis_current = true;

    free_matrix_data(&svd_input);
    free_matrix_data(&slpk_rand_proj);
    free_matrix_data(QRmgr.A);
    free(QRmgr.tau);
    release_context(&slpk_rand_proj);

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
