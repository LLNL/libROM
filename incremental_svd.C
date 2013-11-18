#include "incremental_svd.h"

#include "mpi.h"

#include <string.h>
#include <stdio.h>

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
      d_rank = 0;
      d_size = 1;
   }
}

incremental_svd::~incremental_svd()
{
   // Delete data members.
   for (int i = 0; i < d_num_time_intervals; ++i) {
      if (d_U[i]) {
         delete d_U[i];
      }
      if (d_L[i]) {
         delete d_L[i];
      }
      if (d_S[i]) {
         delete d_S[i];
      }
      if (d_model[i]) {
         delete d_model[i];
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

#ifdef DEBUG_ROMS
   if (d_rank == 0) {
      double* U = new double [d_dim*d_num_increments];
      for (int i = 0; i < d_num_time_intervals; ++i) {
         // Print d_S[i].
         const Matrix* this_d_S = d_S[i];
         for (int row = 0; row < d_num_increments; ++row) {
            for (int col = 0; col < d_num_increments; ++col) {
               printf("%g  ", this_d_S->item(row, col));
            }
            printf("\n");
         }
         printf("\n");

         // Print d_L[i].
         const Matrix* this_d_L = d_L[i];
         for (int row = 0; row < d_num_increments; ++row) {
            for (int col = 0; col < d_num_increments; ++col) {
               printf("%g  ", this_d_L->item(row, col));
            }
            printf("\n");
         }
         printf("\n");

         // Print process 0's part of d_U[i].
         const Matrix* this_d_U = d_U[i];
         for (int row = 0; row < d_dim; ++row) {
            for (int col = 0; col < d_num_increments; ++col) {
               printf("%.16e ", this_d_U->item(row, col));
            }
            printf("\n");
         }

         // Gather other processor's parts of d_U[i] and print them.
         for (int proc = 1; proc < d_size; ++proc) {
            MPI_Status status;
            MPI_Recv(U, d_dim*d_num_increments, MPI_DOUBLE, proc,
                     COMMUNICATE_U, MPI_COMM_WORLD, &status);
            int idx = 0;
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
         MPI_Isend(&d_U[i]->item(0, 0), d_dim*d_num_increments, MPI_DOUBLE, 0,
                   COMMUNICATE_U, MPI_COMM_WORLD, &request);
      }
   }
#endif
}

const Matrix*
incremental_svd::getModel(
   double time)
{
   CAROM_ASSERT(0 < d_num_time_intervals);
   int i;
   for (i = 0; i < d_num_time_intervals-1; ++i) {
      if (d_time_interval_start_times[i] <= time &&
          d_time_interval_start_times[i+1] < time) {
         break;
      }
   }

   // If this model is for the last time interval then it may not be up to
   // date so recompute it.
   if (i == d_num_time_intervals-1) {
      if (d_model[i] != 0) {
         delete d_model[i];
      }
      d_model[i] = d_U[i]->Mult(*d_L[i]);
   }
   else {
      CAROM_ASSERT(d_model[i] != 0);
   }
   return d_model[i];
}

void
incremental_svd::writeModel(
   const std::string& base_file_name)
{
   CAROM_ASSERT(!base_file_name.empty());

   char tmp[10];
   sprintf(tmp, ".%06d", d_rank);
   std::string full_file_name = base_file_name + tmp;
   database.create(full_file_name);
   database.putInteger("num_time_intervals", d_num_time_intervals);
   for (int i = 0; i < d_num_time_intervals; ++i) {
      const Matrix* model = getModel(d_time_interval_start_times[i]);
      database.putDouble("time", d_time_interval_start_times[i]);
      int num_rows = model->numRows();
      database.putInteger("num_rows", num_rows);
      int num_cols = model->numColumns();
      database.putInteger("num_cols", num_cols);
      database.putDoubleArray("model", &model->item(0, 0), num_rows*num_cols);
   }
   database.close();
}

void
incremental_svd::readModel(
   const std::string& base_file_name)
{
   CAROM_ASSERT(!base_file_name.empty());

   char tmp[10];
   sprintf(tmp, ".%06d", d_rank);
   std::string full_file_name = base_file_name + tmp;
   database.open(full_file_name);
   database.getInteger("num_time_intervals", d_num_time_intervals);
   d_time_interval_start_times.resize(d_num_time_intervals);
   d_U.resize(d_num_time_intervals, 0);
   d_L.resize(d_num_time_intervals, 0);
   d_S.resize(d_num_time_intervals, 0);
   d_model.resize(d_num_time_intervals, 0);
   for (int i = 0; i < d_num_time_intervals; ++i) {
      database.getDouble("time", d_time_interval_start_times[i]);
      int num_rows;
      database.getInteger("num_rows", num_rows);
      int num_cols;
      database.getInteger("num_cols", num_cols);
      d_model[i] = new Matrix(num_rows, num_cols, true, d_rank, d_size);
      database.getDoubleArray("model",
                              &d_model[i]->item(0, 0),
                              num_rows*num_cols);
   }
   database.close();
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

   // Build d_S for this new time interval.
   d_S[d_num_time_intervals-1] = new Matrix(1, 1, false, d_rank, d_size);
   Vector u_vec(u, d_dim, true, d_rank, d_size);
   double norm_u = u_vec.norm();
   d_S[d_num_time_intervals-1]->item(0, 0) = norm_u;

   // Build d_L for this new time interval.
   d_L[d_num_time_intervals-1] = new Matrix(1, 1, false, d_rank, d_size);
   d_L[d_num_time_intervals-1]->item(0, 0) = 1.0;

   // Build d_U for this new time interval.
   d_U[d_num_time_intervals-1] = new Matrix(d_dim, 1, true, d_rank, d_size);
   Matrix* this_d_U = d_U[d_num_time_intervals-1];
   for (int i = 0; i < d_dim; ++i) {
      this_d_U->item(i, 0) = u[i]/norm_u;
   }

   // Set the model for this time interval to 0.
   d_model[d_num_time_intervals-1] = 0;

   // If this is not the first time interval then compute the model parameters
   // for the previous time interval and delete the now unnecessary storage for
   // d_U, d_L, and d_S for the previous time interval.
   if (d_num_time_intervals > 1) {
      if (d_model[d_num_time_intervals-2] != 0) {
         delete d_model[d_num_time_intervals-2];
      }
      d_model[d_num_time_intervals-2] =
         d_U[d_num_time_intervals-2]->Mult(*d_L[d_num_time_intervals-2]);
      delete d_U[d_num_time_intervals-2];
      d_U[d_num_time_intervals-2] = 0;
      delete d_L[d_num_time_intervals-2];
      d_L[d_num_time_intervals-2] = 0;
      delete d_S[d_num_time_intervals-2];
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
   Vector* j;
   Matrix* P;
   compute_J_P_normJ(u, j, P);

   // Compute cm = u.u and use it to see if this increment is new.
   Vector u_vec(u, d_dim, true, d_rank, d_size);
   double cm = u_vec.dot(u_vec);
   double eps_squared = d_epsilon*d_epsilon;
   bool is_new_increment = ((cm > eps_squared) &&
      ((d_norm_j*d_norm_j)/(eps_squared+cm)) > eps_squared);

   // If this increment is not new and we are skipping redundant increments
   // then clean up and return.
   double norm_j_for_constructQ = d_norm_j;
   if (!is_new_increment) {
      if (d_skip_redundant) {
         delete j;
         delete P;
         return;
      }
      else {
         norm_j_for_constructQ = 0.0;
      }
   }

   // l = P' * u
   Vector* l = P->TransposeMult(u_vec);

   // Done with P.
   delete P;

   // Create Q.
   double* Q;
   constructQ(Q, l, norm_j_for_constructQ);

   // Done with l.
   delete l;

   // Now get the singular value decomposition of Q.
   Matrix* A;
   Matrix* sigma;
   Matrix* B;
   svd(Q, A, sigma, B);

   // Done with Q.
   delete [] Q;

   // At this point we either have a new or a redundant increment and we are
   // not skipping redundant increments.
   if (!is_new_increment && !d_skip_redundant) {
      // This increment is not new and we are not skipping redundant
      // increments.
      addRedundantIncrement(A, sigma);
      delete sigma;
   }
   else if (is_new_increment) {
      // This increment is new.
      // addNewIncrement will assign sigma to d_S[d_num_time_intervals-1]
      // hence it should not be deleted upon return.
      addNewIncrement(j, A, sigma);
   }

   // Clean up.
   delete j;
   delete A;
   delete B;
}

void
incremental_svd::compute_J_P_normJ(
   const double* u,
   Vector*& j,
   Matrix*& P)
{
   // j = u
   j = new Vector(d_dim, true, d_rank, d_size);
   for (int i = 0; i < d_dim; ++i) {
      j->item(i) = u[i];
   }

   // P = d_U[d_num_time_intervals-1] * d_L[d_num_time_intervals-1]
   P = d_U[d_num_time_intervals-1]->Mult(*d_L[d_num_time_intervals-1]);

   // Use modified Gram-Schmidt orthogonalization to modify j.
   orthogonalizeJAndComputeNorm(j, P);
}

void
incremental_svd::orthogonalizeJAndComputeNorm(
   Vector* j,
   const Matrix* P)
{
   double tmp;
   for (int col = 0; col < d_num_increments; ++col) {
      double factor = 0.0;
      tmp = 0.0;
      for (int i = 0; i < d_dim; ++i) {
         tmp += P->item(i, col)*j->item(i);
      }
      if (d_size > 1) {
         MPI_Allreduce(&tmp, &factor, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
      else {
         factor = tmp;
      }

      for (int i = 0; i < d_dim; ++i) {
         j->item(i) -= factor*P->item(i, col);
      }
   }

   // Normalize j.
   d_norm_j = j->normalize();
}

void
incremental_svd::constructQ(
   double*& Q,
   const Vector* l,
   double norm_j)
{
   // Create Q.
   Q = new double [(d_num_increments+1)*(d_num_increments+1)];

   // Fill Q in column major order.
   int q_idx = 0;
   const Matrix* this_d_S = d_S[d_num_time_intervals-1];
   for (int row = 0; row < d_num_increments; ++row) {
      q_idx = row;
      for (int col = 0; col < d_num_increments; ++col) {
         Q[q_idx] = this_d_S->item(row, col);
         q_idx += d_num_increments+1;
      }
      Q[q_idx] = l->item(row);
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
   Matrix*& U,
   Matrix*& S,
   Matrix*& V)
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
   V = new Matrix(d_num_increments+1,
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
   delete [] work;

   // Place sigma into S.
   for (int i = 0; i < d_num_increments+1; ++i) {
      S->item(i, i) = sigma[i];
   }
   delete [] sigma;

   // U and V are in column major order.  Convert them to row major order.
   for (int row = 0; row < d_num_increments+1; ++row) {
      for (int col = row+1; col < d_num_increments+1; ++col) {
         double tmp = U->item(row, col);
         U->item(row, col) = U->item(col, row);
         U->item(col, row) = tmp;
         tmp = V->item(row, col);
         V->item(row, col) = V->item(col, row);
         V->item(col, row) = tmp;
      }
   }
}

void
incremental_svd::addRedundantIncrement(
   const Matrix* A,
   const Matrix* sigma)
{
   // Chop a row and a column off of A to form Amod.  Also form
   // d_S[d_num_time_intervals-1] by chopping a row and a column off of sigma.
   Matrix Amod(d_num_increments, d_num_increments, false, d_rank, d_size);
   Matrix* this_d_S = d_S[d_num_time_intervals-1];
   for (int row = 0; row < d_num_increments; ++row){
      for (int col = 0; col < d_num_increments; ++col) {
         Amod.item(row, col) = A->item(row, col);
         this_d_S->item(row, col) = sigma->item(row, col);
      }
   }

   // Multiply d_L[d_num_time_intervals-1] and Amod and put result into
   // d_L[d_num_time_intervals-1].
   Matrix* L_times_Amod = d_L[d_num_time_intervals-1]->Mult(Amod);
   delete d_L[d_num_time_intervals-1];
   d_L[d_num_time_intervals-1] = L_times_Amod;
}

void
incremental_svd::addNewIncrement(
   const Vector* j,
   const Matrix* A,
   Matrix* sigma)
{
   // Add j as a new column of d_U[d_num_time_intervals-1].
   const Matrix* this_d_U = d_U[d_num_time_intervals-1];
   Matrix* newU = new Matrix(d_dim, d_num_increments+1, true, d_rank, d_size);
   for (int row = 0; row < d_dim; ++row) {
      for (int col = 0; col < d_num_increments; ++col) {
         newU->item(row, col) = this_d_U->item(row, col);
      }
      newU->item(row, d_num_increments) = j->item(row);
   }
   delete d_U[d_num_time_intervals-1];
   d_U[d_num_time_intervals-1] = newU;

   // Form a temporary matrix, tmp, by add another row and column to
   // d_L[d_num_time_intervals-1].  Only the last value in the new row/column
   // is non-zero and it is 1.
   Matrix tmp(d_num_increments+1, d_num_increments+1, false, d_rank, d_size);
   const Matrix* this_d_L = d_L[d_num_time_intervals-1];
   for (int row = 0; row < d_num_increments; ++row) {
      for (int col = 0; col < d_num_increments; ++col) {
         tmp.item(row, col) = this_d_L->item(row, col);
      }
      tmp.item(row, d_num_increments) = 0.0;
   }
   for (int col = 0; col < d_num_increments; ++col) {
      tmp.item(d_num_increments, col) = 0.0;
   }
   tmp.item(d_num_increments, d_num_increments) = 1.0;

   // d_L[d_num_time_intervals-1] = tmp*A
   delete d_L[d_num_time_intervals-1];
   d_L[d_num_time_intervals-1] = tmp.Mult(*A);

   // d_S[d_num_time_intervals-1] = sigma.
   delete d_S[d_num_time_intervals-1];
   d_S[d_num_time_intervals-1] = sigma;

   // We now have another increment.
   ++d_num_increments;
}

}
