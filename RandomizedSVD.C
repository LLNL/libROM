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

   /*

     %   "Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate
     %    Matrix Decompositions" by N. Halko, P. G. Martinsson, and J. A. Tropp

     [m,n] = size(A);
     if nargin == 1
         k = (m>n)*n + (m<=n)*m;
     end

     if m > n
         Omega = randn(n,k);
         Y = A * Omega;
         [Q,~] = qr(Y,0);
         B = Q'*A;
         [Uhat,Sigma,V] = svd(B);
         U = Q*Uhat;
     else
         Omega = randn(m,k);
         Y = A' * Omega;
         [Q,~] = qr(Y,0);
         B = Q'*A';
         [Uhat,Sigma,Vd] = svd(B);
         Ud = Q*Uhat;
         U = Vd;
         V = Ud;
     end


   */

   void
   RandomizedSVD::computeSVD()
   {
       d_samples->n = d_num_samples;
       delete_factorizer();

       int num_rows = d_total_dim;
       int num_cols = d_num_samples;
       int proj_rows;
       if (d_subspace_dim < 1 || d_subspace_dim > std::min(num_rows, num_cols)) {
         d_subspace_dim = std::min(num_rows, num_cols);
       }
       d_subspace_dim = 8;
       int proj_cols = d_subspace_dim;
       const char* prefix = "";

       // Get snapshot matrix in distributed format.
       // If there are less dimensions than samples, use the transpose instead.
       Matrix* snapshot_matrix;
       if (num_rows > num_cols) {
         snapshot_matrix = new Matrix(d_dims[static_cast<unsigned>(d_rank)],
             d_num_samples, true);
             for (int rank = 0; rank < d_num_procs; ++rank) {
                // gather_transposed_block does the same as gather_block, but transposes
                // it; here, it is used to go from column-major to row-major order.
                gather_transposed_block(&snapshot_matrix->item(0, 0), d_samples.get(),
                                       d_istarts[static_cast<unsigned>(rank)] + 1, 1,
                                       d_dims[static_cast<unsigned>(rank)], d_num_samples,
                                       rank);
             }
        }
        else {
          snapshot_matrix = new Matrix(d_num_samples,
             d_dims[static_cast<unsigned>(d_rank)], true);
             for (int rank = 0; rank < d_num_procs; ++rank) {
                gather_block(&snapshot_matrix->item(0, 0), d_samples.get(),
                                    d_istarts[static_cast<unsigned>(rank)] + 1, 1,
                                    d_dims[static_cast<unsigned>(rank)], d_num_samples,
                                    rank);
             }
        }

       Matrix* rand_mat;
       if (num_rows > num_cols) {
         rand_mat = new Matrix(snapshot_matrix->numColumns(), d_subspace_dim, false, true);
         proj_rows = snapshot_matrix->numColumns();
       }
       else {
         rand_mat = new Matrix(snapshot_matrix->numRows(), d_subspace_dim, false, true);
         proj_rows = snapshot_matrix->numRows();
       }

       // Project snapshot matrix onto random subspace
       Matrix* rand_proj;
       rand_proj = snapshot_matrix->mult(rand_mat);

       // Get QR factorization of random projection
       // Because Fortran is in column-major order, we are passing in the
       // transpose of the snapshot matrix and obtaining the LQ factorization
       // instead.
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
       Matrix* lq_A = new Matrix(d_dims[static_cast<unsigned>(d_rank)],
           d_subspace_dim, true);
       Matrix* svd_input_mat = new Matrix(d_dims[static_cast<unsigned>(d_rank)],
           d_num_samples, true);
       for (int rank = 0; rank < d_num_procs; ++rank) {
         
          gather_transposed_block(&lq_A->item(0, 0), QRmgr.A,
                                 d_istarts[static_cast<unsigned>(rank)] + 1, 1,
                                 d_dims[static_cast<unsigned>(rank)], d_subspace_dim,
                                 rank);
       }

       for (int i = 0; i < lq_A->numColumns(); i++) {
           for (int j = 0; j < i; j++) {
             lq_A->item(j, i) = 0;
           }
           lq_A->item(i, i) = 1;
       }

       Matrix* Q_elementary_household_reflectors = new Matrix(*lq_A);
       int snapshot_matrix_num_cols = snapshot_matrix->numColumns();
       for (int i = 0; i < snapshot_matrix_num_cols; i++) {
         Vector* snapshot_matrix_col = snapshot_matrix->getColumn(i);
         for (int j = Q_elementary_household_reflectors->numColumns() - 1; j >= 0; j--) {
           Vector* v = Q_elementary_household_reflectors->getColumn(j);
           Vector* vvw = NULL;
           v->mult(v->inner_product(snapshot_matrix_col), vvw);
           Vector* tauvvw = NULL;
           vvw->mult(QRmgr.tau[j], tauvvw);
           Vector* temp = new Vector(*snapshot_matrix_col);
           temp->minus(*tauvvw, snapshot_matrix_col);
           delete v, vvw, tauvvw, temp;
         }
         for (int k = 0; k < Q_elementary_household_reflectors->numRows(); k++) {
           svd_input_mat->item(k, i) = snapshot_matrix_col->item(k);
         }
       }

       SLPK_Matrix svd_input;
       initialize_matrix(&svd_input, svd_input_mat->numColumns(), svd_input_mat->numDistributedRows(),
         d_nprow, d_npcol, d_blocksize, d_blocksize);

       for (int rank = 0; rank < d_num_procs; ++rank) {
          scatter_block(&svd_input, 1, row_offset[rank] + 1,
                    svd_input_mat->getData(), svd_input_mat->numColumns(),
                    row_offset[rank + 1] - row_offset[rank], rank);
       }

      // This block does the actual ScaLAPACK call to do the factorization.
      svd_init(d_factorizer.get(), &svd_input);

      d_factorizer->dov = 1;
      factorize(d_factorizer.get());

      // Compute how many basis vectors we will actually return.
      int sigma_cutoff = 0, hard_cutoff = d_num_samples;
      if (d_singular_value_tol == 0) {
         sigma_cutoff = std::numeric_limits<int>::max();
      } else {
         for (int i = 0; i < d_num_samples; ++i) {
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
      d_basis = new Matrix(d_dim, ncolumns, true);
      d_S = new Matrix(ncolumns, ncolumns, false);
      {
         CAROM_VERIFY(ncolumns >= 0);
         unsigned nc = static_cast<unsigned>(ncolumns);
         memset(&d_S->item(0, 0), 0, nc*nc*sizeof(double));
      }

      d_basis_right = new Matrix(ncolumns, d_subspace_dim, false);

      // Since the input to the SVD was transposed, U and V are switched.
      for (int rank = 0; rank < d_num_procs; ++rank) {
         // V is computed in the transposed order so no reordering necessary.
         gather_block(&d_basis->item(0, 0), d_factorizer->V,
                                 1, d_istarts[static_cast<unsigned>(rank)]+1,
                                 ncolumns, d_dims[static_cast<unsigned>(rank)],
                                 rank);

         // gather_transposed_block does the same as gather_block, but transposes
         // it; here, it is used to go from column-major to row-major order.
         gather_transposed_block(&d_basis_right->item(0, 0), d_factorizer->U, 1, 1,
                      ncolumns, d_subspace_dim, rank);
      }
      for (int i = 0; i < ncolumns; ++i)
         d_S->item(i, i) = d_factorizer->S[static_cast<unsigned>(i)];
      d_this_interval_basis_current = true;

      Matrix* new_d_basis = new Matrix(*d_basis);
      for (int i = 0; i < d_basis->numColumns(); i++) {
        Vector* basis_col = d_basis->getColumn(i);
         for (int j = 0; j < Q_elementary_household_reflectors->numColumns(); j++) {
          Vector* v = Q_elementary_household_reflectors->getColumn(j);
          Vector* vvw = NULL;
          v->mult(v->inner_product(basis_col), vvw);
          Vector* tauvvw = NULL;
          vvw->mult(QRmgr.tau[j], tauvvw);
          Vector* temp = new Vector(*basis_col);
          temp->minus(*tauvvw, basis_col);
          delete v, vvw, tauvvw, temp;
        }
        for (int k = 0; k < d_basis->numRows(); k++) {
          new_d_basis->item(k, i) = basis_col->item(k);
        }
        delete basis_col;
      }
      free_matrix_data(&svd_input);
      free_matrix_data(&slpk_rand_proj);
      free_matrix_data(QRmgr.A);
      release_context(&slpk_rand_proj);
      release_context(&svd_input);
      free(QRmgr.tau);
      delete d_basis;
      if (num_rows > num_cols) {
        d_basis = new_d_basis;
      }
      else {
        d_basis = d_basis_right;
        d_basis_right = new_d_basis;
      }

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
