#include "incremental_svd.h"
#include "vector_utils.h"
#include <assert.h>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <mpi.h>

extern "C" {
void dgesdd_(char*, int*, int*, double*, int*,
             double*, double*, int*, double*, int*,
             double*, int*, int*, int*);
}

// #define DEBUG

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
   d_U(0),
   d_L(0),
   d_S(0),
   d_model(0),
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
      d_rank = 1;
      d_size = 1;
   }
}

incremental_svd::~incremental_svd()
{
   // Delete data members.
   for (int i = 0; i < d_num_time_intervals; ++i) {
      if (d_U[i]) {
         delete [] d_U[i];
      }
      if (d_L[i]) {
         delete [] d_L[i];
      }
      if (d_S[i]) {
         delete [] d_S[i];
      }
      if (d_model[i]) {
         delete [] d_model[i];
      }
   }
}

void
incremental_svd::increment(
   double* u_in,
   double time)
{
   // If this is the first SVD then build it.  Otherwise add this increment to
   // the system.
   if (d_num_increments == 0 ||
       d_num_increments >= d_increments_per_time_interval) {
     buildInitialSVD(u_in, time);
   }
   else {
      buildIncrementalSVD(u_in);
   }

#ifdef DEBUG
   if (d_rank == 0) {
      double* U = new double [d_dim*d_num_increments];
      for (int i = 0; i < d_num_time_intervals; ++i) {
         // Print d_S[i].
         const double* this_d_S = d_S[i];
         int idx = 0;
         for (int row = 0; row < d_num_increments; ++row) {
            for (int col = 0; col < d_num_increments; ++col) {
               printf("%g  ", this_d_S[idx++]);
            }
            printf("\n");
         }
         printf("\n");

         // Print d_L[i].
         const double* this_d_L = d_L[i];
         idx = 0;
         for (int row = 0; row < d_num_increments; ++row) {
            for (int col = 0; col < d_num_increments; ++col) {
               printf("%g  ", this_d_L[idx++]);
            }
            printf("\n");
         }
         printf("\n");

         // Print process 0's part of d_U[i].
         const double* this_d_U = d_U[i];
         idx = 0;
         for (int row = 0; row < d_dim; ++row) {
            for (int col = 0; col < d_num_increments; ++col) {
               printf("%.16e ", this_d_U[idx++]);
            }
            printf("\n");
         }

         // Gather other processor's parts of d_U[i] and print them.
         for (int proc = 1; proc < d_size; ++proc) {
            MPI_Status status;
            MPI_Recv(U, d_dim*d_num_increments, MPI_DOUBLE, proc,
                     COMMUNICATE_U, MPI_COMM_WORLD, &status);
            idx = 0;
            for (int row = 0; row < d_dim; ++row) {
               for (int col = 0; col < d_num_increments; ++col) {
                  printf("%.16e ", U[idx++]);
               }
               printf("\n");
            }
         }
         printf("============================================================\n");
      }
      delete [] U;
   }
   else {
      // DOES THIS WORK??????????????????????????????????????
      // Send this processor's part of d_U[i] to process 0.
      for (int i = 0; i < d_num_time_intervals; ++i) {
         MPI_Request request;
         MPI_Isend(d_U[i], d_dim*d_num_increments, MPI_DOUBLE, 0,
                   COMMUNICATE_U, MPI_COMM_WORLD, &request);
      }
   }
#endif
}

double*
incremental_svd::getModel(
   double time) const
{
   assert(0 < d_num_time_intervals);
   int i;
   for (i = 0; i < d_num_time_intervals-1; ++i) {
      if (d_time_interval_start_times[i] <= time &&
          d_time_interval_start_times[i+1] < time) {
         break;
      }
   }
   if (i == d_num_time_intervals-1) {
      return DistributedMatLocalMatMult(d_U[i], d_dim, d_num_increments,
                                        d_L[i], d_num_increments,
                                        d_num_increments);
   }
   else {
      assert(d_model[i] != 0);
      return d_model[i];
   }
}

void
incremental_svd::buildInitialSVD(
   double* u,
   double time)
{
   // We have a new time interval.
   ++d_num_time_intervals;
   d_time_interval_start_times.resize(d_num_time_intervals);
   d_time_interval_start_times[d_num_time_intervals-1] = time;
   d_U.resize(d_num_time_intervals);
   d_L.resize(d_num_time_intervals);
   d_S.resize(d_num_time_intervals);
   d_model.resize(d_num_time_intervals);

   // Build d_S[d_num_time_intervals-1].
   d_S[d_num_time_intervals-1] = new double [1];
   double norm_u = norm(u, d_dim, d_size);
   d_S[d_num_time_intervals-1][0] = norm_u;

   // Build d_L[d_num_time_intervals-1].
   d_L[d_num_time_intervals-1] = new double [1];
   d_L[d_num_time_intervals-1][0] = 1.0;

   // Build d_U[d_num_time_intervals-1].
   d_U[d_num_time_intervals-1] = new double [d_dim];
   double* this_d_U = d_U[d_num_time_intervals-1];
   for (int i = 0; i < d_dim; ++i) {
      this_d_U[i] = u[i]/norm_u;
   }

   // If this is not the first time interval then compute the model parameters
   // for the previous time interval and delete the now unnecessary storage for
   // d_U, d_L, and d_S for the previous time interval.
   if (d_num_time_intervals > 1) {
      if (d_model[d_num_time_intervals-2] != 0) {
         delete [] d_model[d_num_time_intervals-2];
      }
      d_model[d_num_time_intervals-2] =
         DistributedMatLocalMatMult(d_U[d_num_time_intervals-2],
                                    d_dim,
                                    d_num_increments,
                                    d_L[d_num_time_intervals-2],
                                    d_num_increments,
                                    d_num_increments);
      delete [] d_U[d_num_time_intervals-2];
      d_U[d_num_time_intervals-2] = 0;
      delete [] d_L[d_num_time_intervals-2];
      d_L[d_num_time_intervals-2] = 0;
      delete [] d_S[d_num_time_intervals-2];
      d_S[d_num_time_intervals-2] = 0;
   }

   // We now have the first increment for the new time interval.
   d_num_increments = 1;
}

void
incremental_svd::buildIncrementalSVD(
   double* u)
{
   // Compute j, P and the norm of j.
   double* j;
   double* P;
   compute_J_P_normJ(u, j, P);

   // l = P' * u
   double* l =
      DistributedMatTransposeDistributedMatMult(P, d_dim, d_num_increments,
                                                u, d_dim, 1, d_size);

   // Compute cm = u.u and use it to see if this increment is new.
   double cm = inner_product(u, u, d_dim, d_size);
   double eps_squared = d_epsilon*d_epsilon;
   bool is_new_increment = ((cm > eps_squared) &&
      ((d_norm_j*d_norm_j)/(eps_squared+cm)) > eps_squared);

   // If this increment is not new and we are skipping redundant increments
   // then clean up and return.
   double norm_j_for_constructQ = d_norm_j;
   if (!is_new_increment) {
      if (d_skip_redundant) {
         delete [] j;
         delete [] l;
         return;
      }
      else {
         norm_j_for_constructQ = 0.0;
      }
   }

   // Create Q.
   double* Q;
   constructQ(Q, l, norm_j_for_constructQ);

   // Done with l.
   delete [] l;

   // Now get the singular value decomposition of Q.
   double* A;
   double* sigma;
   double* B;
   svd(Q, A, sigma, B);

   // Done with Q.
   delete [] Q;

   // At this point we either have a new or a redundant increment and we are
   // not skipping redundant increments.
   if (!is_new_increment && !d_skip_redundant) {
      // This increment is not new and we are not skipping redundant
      // increments.
      addRedundantIncrement(A, sigma);
      delete [] sigma;
   }
   else if (is_new_increment) {
      // This increment is new.
      // addNewIncrement will assign sigma to d_S[d_num_time_intervals-1]
      // hence it should not be deleted upon return.
      addNewIncrement(j, A, sigma);
   }

   // Clean up.
   delete [] P;
   delete [] j;
   delete [] A;
   delete [] B;
}

void
incremental_svd::compute_J_P_normJ(
   double* u,
   double*& j,
   double*& P)
{
   // j = u
   j = new double [d_dim];
   memcpy(j, u, d_dim*sizeof(double));

   // P = d_U[d_num_time_intervals-1] * d_L[d_num_time_intervals-1]
   P = DistributedMatLocalMatMult(d_U[d_num_time_intervals-1],
                                  d_dim,
                                  d_num_increments,
                                  d_L[d_num_time_intervals-1],
                                  d_num_increments,
                                  d_num_increments);

   // Use modified Gram-Schmidt orthogonalization to modify j.
   orthogonalizeJAndComputeNorm(j, P);
}

void
incremental_svd::orthogonalizeJAndComputeNorm(
   double* j,
   double* P)
{
   double tmp;
   for (int col = 0; col < d_num_increments; ++col) {
      double factor = 0.0;
      tmp = 0.0;
      int Pidx = col;
      for (int i = 0; i < d_dim; ++i) {
         tmp += P[Pidx]*j[i];
         Pidx += d_num_increments;
      }
      if (d_size > 1) {
         MPI_Allreduce(&tmp, &factor, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
      else {
         factor = tmp;
      }
      
      Pidx = col;
      for (int i = 0; i < d_dim; ++i) {
         j[i] -= factor*P[Pidx];
         Pidx += d_num_increments;
      }
   }

   // Normalize j.
   d_norm_j = normalize(j, d_dim, d_size);
}

void
incremental_svd::constructQ(
   double*& Q,
   double* l,
   double norm_j)
{
   // Create Q.
   Q = new double [(d_num_increments+1)*(d_num_increments+1)];

   // Fill Q in column major order.
   int q_idx = 0;
   int s_idx = 0;
   const double* this_d_S = d_S[d_num_time_intervals-1];
   for (int row = 0; row < d_num_increments; ++row) {
      q_idx = row;
      for (int col = 0; col < d_num_increments; ++col) {
         Q[q_idx] = this_d_S[s_idx++];
         q_idx += d_num_increments+1;
      }
      Q[q_idx] = l[row];
   }
   q_idx = d_num_increments;
   for (int col = 0; col < d_num_increments; ++col) {
      Q[q_idx] = 0.0;
      q_idx += d_num_increments+1;
   }
   Q[q_idx] = norm_j;
}

void
incremental_svd::svd(
   double* A,
   double*& U,
   double*& S,
   double*& V)
{
   // Construct U, S, and V.
   int mat_size = (d_num_increments+1)*(d_num_increments+1);
   U = new double [mat_size];
   S = new double [mat_size];
   V = new double [mat_size];
   for (int i = 0; i < mat_size; ++i) {
      S[i] = 0.0;
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
           sigma, U, &ldu, V, &ldv,
           work, &lwork, iwork, &info);
   delete [] work;

   // Place sigma into S.
   for (int i = 0; i < d_num_increments+1; ++i) {
      S[i*(d_num_increments+1)+i] = sigma[i];
   }
   delete [] sigma;

   // U and V are in column major order.  Convert them to row major order.
   for (int row = 0; row < d_num_increments+1; ++row) {
      for (int col = row+1; col < d_num_increments+1; ++col) {
         double tmp = U[row*(d_num_increments+1)+col];
         U[row*(d_num_increments+1)+col] = U[col*(d_num_increments+1)+row];
         U[col*(d_num_increments+1)+row] = tmp;
         tmp = V[row*(d_num_increments+1)+col];
         V[row*(d_num_increments+1)+col] = V[col*(d_num_increments+1)+row];
         V[col*(d_num_increments+1)+row] = tmp;
      }
   }
}

void
incremental_svd::addRedundantIncrement(
   double* A,
   double* sigma)
{
   // Chop a row and a column off of A to form Amod.  Also form
   // d_S[d_num_time_intervals-1] by chopping a row and a column off of sigma.
   double* Amod = new double [d_num_increments*d_num_increments];
   int lhs_idx = 0;
   int rhs_idx = 0;
   double* this_d_S = d_S[d_num_time_intervals-1];
   for (int row = 0; row < d_num_increments; ++row){
      for (int col = 0; col < d_num_increments; ++col) {
         Amod[lhs_idx] = A[rhs_idx];
         this_d_S[lhs_idx] = sigma[rhs_idx];
         ++lhs_idx;
         ++rhs_idx;
      }
      ++rhs_idx;
   }

   // Multiply d_L[d_num_time_intervals-1] and Amod and put result into
   // d_L[d_num_time_intervals-1].
   double* L_times_Amod =
      LocalMatLocalMatMult(d_L[d_num_time_intervals-1],
                           d_num_increments,
                           d_num_increments,
                           Amod,
                           d_num_increments,
                           d_num_increments);
   delete [] d_L[d_num_time_intervals-1];
   d_L[d_num_time_intervals-1] = L_times_Amod;

   // Clean up.
   delete [] Amod;
}

void
incremental_svd::addNewIncrement(
   double* j,
   double* A,
   double* sigma)
{
   // Add j as a new column of d_U[d_num_time_intervals-1].
   const double* this_d_U = d_U[d_num_time_intervals-1];
   double* newU = new double [d_dim*(d_num_increments+1)];
   int lhs_idx = 0;
   int rhs_idx = 0;
   for (int row = 0; row < d_dim; ++row) {
      for (int col = 0; col < d_num_increments; ++col) {
         newU[lhs_idx++] = this_d_U[rhs_idx++];
      }
      newU[lhs_idx++] = j[row];
   }
   delete [] d_U[d_num_time_intervals-1];
   d_U[d_num_time_intervals-1] = newU;

   // Add another row and column to d_L[d_num_time_intervals-1].  Only the last
   // value in the new row/column is non-zero and it is 1.
   double* newL = new double [(d_num_increments+1)*(d_num_increments+1)];
   lhs_idx = 0;
   rhs_idx = 0;
   const double* this_d_L = d_L[d_num_time_intervals-1];
   for (int row = 0; row < d_num_increments; ++row) {
      for (int col = 0; col < d_num_increments; ++col) {
         newL[lhs_idx++] = this_d_L[rhs_idx++];
      }
      newL[lhs_idx++] = 0.0;
   }
   for (int col = 0; col < d_num_increments; ++col) {
      newL[lhs_idx++] = 0.0;
   }
   newL[lhs_idx] = 1.0;
   delete [] d_L[d_num_time_intervals-1];
   d_L[d_num_time_intervals-1] =
      LocalMatLocalMatMult(newL, d_num_increments+1, d_num_increments+1,
                           A, d_num_increments+1, d_num_increments+1);
   delete [] newL;

   // d_S[d_num_time_intervals-1] = sigma.
   delete [] d_S[d_num_time_intervals-1];
   d_S[d_num_time_intervals-1] = sigma;

   // We now have another increment.
   ++d_num_increments;
}

}
