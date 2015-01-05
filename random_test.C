/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: A test of both the static and incremental fast update
 *              algorithms and samplers.  Random numbers are generated as the
 *              state vectors in such a way that the global state vector is the
 *              same for all processor decompositions when
 *              dim * number of processors is constant.
 *
 *****************************************************************************/

#include "IncrementalSVDBasisGenerator.h"
#include "StaticSVDBasisGenerator.h"

#include "mpi.h"

#include <stdio.h>

CAROM::Matrix*
transposeMult(
   const CAROM::Matrix* umat,
   const CAROM::Matrix* dmat)
{
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   int umat_num_cols = umat->numColumns();
   int dmat_num_cols = dmat->numColumns();
   int dmat_num_rows = dmat->numRows();
   int new_mat_size = umat_num_cols*dmat_num_cols;
   CAROM::Matrix* local_result = new CAROM::Matrix(umat_num_cols,
      dmat_num_cols,
      false);
   for (int umat_col = 0; umat_col < umat_num_cols; ++umat_col) {
      for (int dmat_col = 0; dmat_col < dmat_num_cols; ++dmat_col) {
         int umat_row = dmat_num_rows*rank;
         local_result->item(umat_col, dmat_col) = 0.0;
         for (int entry = 0; entry < dmat_num_rows; ++entry) {
            local_result->item(umat_col, dmat_col) +=
               umat->item(umat_row, umat_col)*dmat->item(entry, dmat_col);
            ++umat_row;
         }
      }
   }
   CAROM::Matrix* result;
   int num_procs;
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   if (num_procs > 1) {
      result = new CAROM::Matrix(umat_num_cols, dmat_num_cols, false);
      MPI_Allreduce(&local_result->item(0, 0),
         &result->item(0, 0),
         new_mat_size,
         MPI_DOUBLE,
         MPI_SUM,
         MPI_COMM_WORLD);
      delete local_result;
   }
   else {
      result = local_result;
   }
   return result;
}

int main(int argc, char* argv[])
{
   if (argc != 4) {
      printf("Usage:  random_test dim num_snapshots num_lin_dep_snapshots\n");
      return 1;
   }
   int dim = atoi(argv[1]);
   int num_snapshots = atoi(argv[2]);
   int num_lin_dep_snapshots = atoi(argv[3]);
   int num_lin_indep_snapshots = num_snapshots - num_lin_dep_snapshots;
   MPI_Init(&argc, &argv);
   CAROM::IncrementalSVDBasisGenerator inc_basis_generator(dim,
      1.0e-6,
      false,
      num_snapshots,
      1.0e-2,
      0.001,
      true,
      "");
   CAROM::StaticSVDBasisGenerator static_basis_generator(dim,
      num_snapshots,
      "");
   int size;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   srand(1);
   for (int i = 0; i < dim*num_snapshots*rank; ++i) {
      double random = rand();
      random = random/RAND_MAX;
   }
   double** M = new double* [num_snapshots];
   for (int i = 0; i < num_snapshots; ++i) {
      M[i] = new double [dim];
   }
   for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < num_snapshots; ++j) {
         double random = rand();
         random = random/RAND_MAX;
         M[j][i] = random;
      }
   }
   for (int i = 0; i < dim*num_snapshots*(size-rank-1); ++i) {
      double random = rand();
      random = random/RAND_MAX;
   }
   for (int i = 0; i < num_lin_dep_snapshots; ++i) {
      int col = num_snapshots - i - 1;
      for (int j = 0; j < dim; ++j) {
         M[col][j] = 0;
      }
      for (int j = 0; j < num_lin_indep_snapshots; ++j) {
         double random = rand();
         random = random/RAND_MAX;
         for (int k = 0; k < dim; ++k) {
            M[col][k] += random*M[j][k];
         }
      }
   }
   for (int i = 0; i < num_snapshots; ++i) {
      if (inc_basis_generator.isNextSnapshot(0.01*i)) {
         inc_basis_generator.takeSnapshot(M[i], 0.01*i);
         inc_basis_generator.computeNextSnapshotTime(M[i], M[i], 0.01*i);
      }
   }
   inc_basis_generator.endSnapshots();
   const CAROM::Matrix* inc_basis = inc_basis_generator.getBasis();
   for (int i = 0; i < num_snapshots; ++i) {
      if (static_basis_generator.isNextSnapshot(0.01*i)) {
         static_basis_generator.takeSnapshot(M[i], 0.01*i);
         static_basis_generator.computeNextSnapshotTime(M[i], M[i], 0.01*i);
      }
   }
   static_basis_generator.endSnapshots();
   const CAROM::Matrix* static_basis = static_basis_generator.getBasis();

   // Compute the product of the tranpose of the static basis and the
   // incremental basis.  This should be a unitary matrix.
   CAROM::Matrix* test = transposeMult(static_basis, inc_basis);
   if (rank == 0) {
      for (int row = 0; row < num_snapshots; ++row) {
         for (int col = 0; col < num_lin_indep_snapshots; ++col) {
            printf("%.16e ", test->item(row, col));
         }
         printf("\n");
      }
   }
   delete test;
   for (int i = 0; i < num_snapshots; ++i) {
      delete [] M[i];
   }
   delete [] M;
   MPI_Finalize();
   return 0;
}
