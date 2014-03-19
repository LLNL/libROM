#include "incremental_svd.h"

#include "mpi.h"

#include <string.h>

extern "C" {
void dgesdd_(char*, int*, int*, double*, int*,
             double*, double*, int*, double*, int*,
             double*, int*, int*, int*);
}

namespace CAROM {

const int incremental_svd::COMMUNICATE_U = 666;

incremental_svd::incremental_svd(
   int dim,
   double epsilon,
   bool skip_redundant,
   int increments_per_time_interval) :
   d_dim(dim),
   d_num_increments(0),
   d_epsilon(epsilon),
   d_skip_redundant(skip_redundant),
   d_increments_per_time_interval(increments_per_time_interval),
   d_S(0),
   d_basis(0),
   d_num_time_intervals(0),
   d_time_interval_start_times(0),
   d_norm_j(0.0)
{
   // Get the rank of this process, and get the number of processors.
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init) {
      MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
      MPI_Comm_size(MPI_COMM_WORLD, &d_size);
   }
   else {
      d_rank = 0;
      d_size = 1;
   }
}

incremental_svd::~incremental_svd()
{
   // Delete data members.
   if (d_S) {
      delete d_S;
   }
   for (int i = 0; i < d_num_time_intervals; ++i) {
      if (d_basis[i]) {
         delete d_basis[i];
      }
   }
}

void
incremental_svd::constructQ(
   double*& Q,
   const Vector* l,
   double k)
{
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
incremental_svd::svd(
   double* A,
   Matrix*& U,
   Matrix*& S)
{
   // Construct U, S, and V.
   U = new Matrix(d_num_increments+1,
                  d_num_increments+1,
                  false,
                  d_rank,
                  d_size);
   S = new Matrix(d_num_increments+1,
                  d_num_increments+1,
                  false,
                  d_rank,
                  d_size);
   Matrix* V = new Matrix(d_num_increments+1,
                          d_num_increments+1,
                          false,
                          d_rank,
                          d_size);
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
   dgesdd_(&jobz, &m, &n, A, &lda,
           sigma, &U->item(0, 0), &ldu, &V->item(0, 0), &ldv,
           work, &lwork, iwork, &info);
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
