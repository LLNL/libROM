/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The abstract incremental SVD algorithm defines algorithm
 *              interface.
 *
 *****************************************************************************/

#include "IncrementalSVD.h"

#include "mpi.h"

extern "C" {
void dgesdd_(char*, int*, int*, double*, int*,
             double*, double*, int*, double*, int*,
             double*, int*, int*, int*);
}

namespace CAROM {

const int IncrementalSVD::COMMUNICATE_U = 666;

IncrementalSVD::IncrementalSVD(
   int dim,
   double redundancy_tol,
   bool skip_redundant,
   int increments_per_time_interval) :
   d_dim(dim),
   d_num_increments(0),
   d_redundancy_tol(redundancy_tol),
   d_skip_redundant(skip_redundant),
   d_increments_per_time_interval(increments_per_time_interval),
   d_S(0),
   d_basis(0),
   d_time_interval_start_times(0),
   d_total_dim(0)
{
   CAROM_ASSERT(dim > 0);
   CAROM_ASSERT(redundancy_tol > 0.0);
   CAROM_ASSERT(increments_per_time_interval > 0);

   // Get the number of processors, the dimensions for each process, and the
   // total dimension.
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init) {
      MPI_Comm_size(MPI_COMM_WORLD, &d_size);
   }
   else {
      d_size = 1;
   }
   d_proc_dims.reserve(d_size);
   if (mpi_init) {
      MPI_Allgather(&d_dim,
         1,
         MPI_INT,
         &d_proc_dims[0],
         1,
         MPI_INT,
         MPI_COMM_WORLD);
   }
   else {
      d_proc_dims[0] = d_dim;
   }
   for (int i = 0; i < d_size; ++i) {
      d_total_dim += d_proc_dims[i];
   }
}

IncrementalSVD::~IncrementalSVD()
{
   // Delete data members.
   if (d_S) {
      delete d_S;
   }
   if (d_basis) {
      delete d_basis;
   }
}

void
IncrementalSVD::constructQ(
   double*& Q,
   const Vector* l,
   double k)
{
   CAROM_ASSERT(l != 0);
   CAROM_ASSERT(l->dim() == d_num_increments);

   // Create Q.
   Q = new double [(d_num_increments+1)*(d_num_increments+1)];

   // Fill Q in column major order.
   int q_idx = 0;
   for (int row = 0; row < d_num_increments; ++row) {
      q_idx = row;
      for (int col = 0; col < d_num_increments; ++col) {
         Q[q_idx] = d_S->item(row, col);
         q_idx += d_num_increments+1;
      }
      Q[q_idx] = l->item(row);
   }
   q_idx = d_num_increments;
   for (int col = 0; col < d_num_increments; ++col) {
      Q[q_idx] = 0.0;
      q_idx += d_num_increments+1;
   }
   Q[q_idx] = k;
}

void
IncrementalSVD::svd(
   double* A,
   Matrix*& U,
   Matrix*& S)
{
   CAROM_ASSERT(A != 0);

   // Construct U, S, and V.
   U = new Matrix(d_num_increments+1, d_num_increments+1, false);
   S = new Matrix(d_num_increments+1, d_num_increments+1, false);
   Matrix* V = new Matrix(d_num_increments+1, d_num_increments+1, false);
   for (int row = 0; row < d_num_increments+1; ++row) {
      for (int col = 0; col < d_num_increments+1; ++col) {
         S->item(row, col) = 0.0;
      }
   }

   // Use lapack's dgesdd_ Fortran function to perform the svd.  As this is
   // Fortran A and all the computed matrices are in column major order.
   double* sigma = new double [d_num_increments+1];
   char jobz = 'A';
   int m = d_num_increments+1;
   int n = d_num_increments+1;
   int lda = d_num_increments+1;
   int ldu = d_num_increments+1;
   int ldv = d_num_increments+1;
   int lwork = m*(4*m + 7);
   double* work = new double [lwork];
   int iwork[8*m];
   int info;
   dgesdd_(&jobz,
      &m,
      &n,
      A,
      &lda,
      sigma,
      &U->item(0, 0),
      &ldu,
      &V->item(0, 0),
      &ldv,
      work,
      &lwork,
      iwork,
      &info);
   CAROM_ASSERT(info == 0);
   delete [] work;

   // Place sigma into S.
   for (int i = 0; i < d_num_increments+1; ++i) {
      S->item(i, i) = sigma[i];
   }
   delete [] sigma;

   // U is column major order so convert it to row major order.
   for (int row = 0; row < d_num_increments+1; ++row) {
      for (int col = row+1; col < d_num_increments+1; ++col) {
         double tmp = U->item(row, col);
         U->item(row, col) = U->item(col, row);
         U->item(col, row) = tmp;
      }
   }
   delete V;
}

}
