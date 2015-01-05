/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The concrete implementation of the incremental SVD algorithm
 *              using Matthew Brand's "fast update" method.
 *
 *****************************************************************************/

#include "IncrementalSVDFastUpdate.h"

#include "mpi.h"

#include <cmath>
#include <limits>
#include <stdio.h>

namespace CAROM {

IncrementalSVDFastUpdate::IncrementalSVDFastUpdate(
   int dim,
   double redundancy_tol,
   bool skip_redundant,
   int increments_per_time_interval,
   bool debug_rom) :
   IncrementalSVD(dim,
      redundancy_tol,
      skip_redundant,
      increments_per_time_interval,
      debug_rom),
   d_U(0),
   d_Up(0)
{
}

IncrementalSVDFastUpdate::~IncrementalSVDFastUpdate()
{
   // Delete data members.
   if (d_U) {
      delete d_U;
   }
   if (d_Up) {
      delete d_Up;
   }
}

void
IncrementalSVDFastUpdate::increment(
   const double* u_in,
   double time)
{
   CAROM_ASSERT(u_in != 0);
   CAROM_ASSERT(time >= 0.0);

   // If this is the first SVD then build it.  Otherwise add this increment to
   // the system.
   if (isNewTimeInterval()) {
      buildInitialSVD(u_in, time);
   }
   else {
      buildIncrementalSVD(u_in);
   }

   if (d_debug_rom) {
      int mpi_init;
      MPI_Initialized(&mpi_init);
      int rank;
      if (mpi_init) {
         MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      }
      else {
         rank = 0;
      }
      const Matrix* basis = getBasis();
      if (rank == 0) {
         // Print d_S.
         for (int row = 0; row < d_num_increments; ++row) {
            for (int col = 0; col < d_num_increments; ++col) {
               printf("%.16e  ", d_S->item(row, col));
            }
            printf("\n");
         }
         printf("\n");

         // Print process 0's part of the basis.
         for (int row = 0; row < d_dim; ++row) {
            for (int col = 0; col < d_num_increments; ++col) {
               printf("%.16e ", basis->item(row, col));
            }
            printf("\n");
         }

         // Gather other processor's parts of the basis and print them.
         for (int proc = 1; proc < d_size; ++proc) {
            double* m = new double[d_proc_dims[proc]*d_num_increments];
            MPI_Status status;
            MPI_Recv(m,
               d_proc_dims[proc]*d_num_increments,
               MPI_DOUBLE,
               proc,
               COMMUNICATE_U,
               MPI_COMM_WORLD,
               &status);
            int idx = 0;
            for (int row = 0; row < d_proc_dims[proc]; ++row) {
               for (int col = 0; col < d_num_increments; ++col) {
                  printf("%.16e ", m[idx++]);
               }
               printf("\n");
            }
            delete [] m;
         }
         printf("============================================================\n");
      }
      else {
         // Send this processor's part of the basis to process 0.
         MPI_Request request;
         MPI_Isend(const_cast<double*>(&basis->item(0, 0)),
            d_dim*d_num_increments,
            MPI_DOUBLE,
            0,
            COMMUNICATE_U,
            MPI_COMM_WORLD,
            &request);
      }
   }
}

const Matrix*
IncrementalSVDFastUpdate::getBasis()
{
   CAROM_ASSERT(d_basis != 0);
   return d_basis;
}

void
IncrementalSVDFastUpdate::buildInitialSVD(
   const double* u,
   double time)
{
   CAROM_ASSERT(u != 0);
   CAROM_ASSERT(time >= 0.0);

   // We have a new time interval.

   // If this is not the first time interval then delete d_basis, d_U, d_Up,
   // and d_S of the just completed time interval.
   int num_time_intervals =
      static_cast<int>(d_time_interval_start_times.size());
   if (num_time_intervals > 0) {
      delete d_basis;
      delete d_U;
      delete d_Up;
      delete d_S;
   }
   d_time_interval_start_times.resize(num_time_intervals+1);
   d_time_interval_start_times[num_time_intervals] = time;

   // Build d_S for this new time interval.
   d_S = new Matrix(1, 1, false);
   Vector u_vec(u, d_dim, true);
   double norm_u = u_vec.norm();
   d_S->item(0, 0) = norm_u;

   // Build d_Up for this new time interval.
   d_Up = new Matrix(1, 1, false);
   d_Up->item(0, 0) = 1.0;

   // Build d_U for this new time interval.
   d_U = new Matrix(d_dim, 1, true);
   for (int i = 0; i < d_dim; ++i) {
      d_U->item(i, 0) = u[i]/norm_u;
   }

   // Set the basis for this time interval.
   d_basis = new Matrix(*d_U);

   // We now have the first increment for the new time interval.
   d_num_increments = 1;
}

void
IncrementalSVDFastUpdate::buildIncrementalSVD(
   const double* u)
{
   CAROM_ASSERT(u != 0);

   // l = basis' * u
   Vector u_vec(u, d_dim, true);
   Vector* l = d_basis->transposeMult(u_vec);

   // basisl = basis * l
   Vector* basisl = d_basis->mult(l);

   // Compute k = sqrt(u.u - 2.0*l.l + basisl.basisl) which is ||u - basisl||.
   double k = u_vec.inner_product(u_vec) - 2.0*l->inner_product(l) +
      basisl->inner_product(basisl);
   if (k <= 0) {
      k = 0;
   }
   else {
      k = sqrt(k);
   }

   // Use k to see if this increment is new.
   double default_tol;
   if (d_num_increments == 1) {
      default_tol = std::numeric_limits<double>::epsilon();
   }
   else {
      default_tol = d_num_increments*std::numeric_limits<double>::epsilon()*d_S->item(0, 0);
   }
   if (d_redundancy_tol < default_tol) {
      d_redundancy_tol = default_tol;
   }
   bool is_new_increment;
   if (k < d_redundancy_tol) {
      k = 0;
      is_new_increment = false;
   }
   else {
      is_new_increment = true;
   }

   // Create Q.
   double* Q;
   constructQ(Q, l, k);

   // Now get the singular value decomposition of Q.
   Matrix* A;
   Matrix* sigma;
   svd(Q, A, sigma);

   // Done with Q.
   delete [] Q;

   // At this point we either have a new or a redundant increment and we are
   // not skipping redundant increments.
   if (!is_new_increment && !d_skip_redundant) {
      // This increment is not new and we are not skipping redundant
      // increments.
      addRedundantIncrement(A, sigma);

      // Reorthogonalize if necessary.
      if (fabs(checkUpOrthogonality()) >
          std::numeric_limits<double>::epsilon()*d_num_increments) {
         reOrthogonalize(d_Up);
      }

      delete sigma;
   }
   else if (is_new_increment) {
      // This increment is new.

      // Compute j
      Vector* j = u_vec.minus(basisl);
      for (int i = 0; i < d_dim; ++i) {
         j->item(i) /= k;
      }

      // addNewIncrement will assign sigma to d_S hence it should not be
      // deleted upon return.
      addNewIncrement(j, A, sigma);
      delete j;

      // Reorthogonalize if necessary.
      long int max_U_dim;
      if (d_total_dim > d_num_increments) {
         max_U_dim = d_total_dim;
      }
      else {
         max_U_dim = d_num_increments;
      }
      if (fabs(checkUOrthogonality()) >
          std::numeric_limits<double>::epsilon()*max_U_dim) {
         reOrthogonalize(d_U);
      }
      if (fabs(checkUpOrthogonality()) >
          std::numeric_limits<double>::epsilon()*d_num_increments) {
         reOrthogonalize(d_Up);
      }
   }
   delete basisl;

   // Compute the basis vectors.
   delete d_basis;
   d_basis = d_U->mult(d_Up);

   // Clean up.
   delete l;
   delete A;
}

void
IncrementalSVDFastUpdate::addRedundantIncrement(
   const Matrix* A,
   const Matrix* sigma)
{
   CAROM_ASSERT(A != 0);
   CAROM_ASSERT(sigma != 0);

   // Chop a row and a column off of A to form Amod.  Also form
   // d_S by chopping a row and a column off of sigma.
   Matrix Amod(d_num_increments, d_num_increments, false);
   for (int row = 0; row < d_num_increments; ++row){
      for (int col = 0; col < d_num_increments; ++col) {
         Amod.item(row, col) = A->item(row, col);
         d_S->item(row, col) = sigma->item(row, col);
      }
   }

   // Multiply d_Up and Amod and put result into d_Up.
   Matrix* Up_times_Amod = d_Up->mult(Amod);
   delete d_Up;
   d_Up = Up_times_Amod;
}

void
IncrementalSVDFastUpdate::addNewIncrement(
   const Vector* j,
   const Matrix* A,
   Matrix* sigma)
{
   CAROM_ASSERT(j != 0);
   CAROM_ASSERT(A != 0);
   CAROM_ASSERT(sigma != 0);

   // Add j as a new column of d_U.
   Matrix* newU = new Matrix(d_dim, d_num_increments+1, true);
   for (int row = 0; row < d_dim; ++row) {
      for (int col = 0; col < d_num_increments; ++col) {
         newU->item(row, col) = d_U->item(row, col);
      }
      newU->item(row, d_num_increments) = j->item(row);
   }
   delete d_U;
   d_U = newU;

   // Form a temporary matrix, tmp, by add another row and column to
   // d_Up.  Only the last value in the new row/column is non-zero and it is 1.
   Matrix tmp(d_num_increments+1, d_num_increments+1, false);
   for (int row = 0; row < d_num_increments; ++row) {
      for (int col = 0; col < d_num_increments; ++col) {
         tmp.item(row, col) = d_Up->item(row, col);
      }
      tmp.item(row, d_num_increments) = 0.0;
   }
   for (int col = 0; col < d_num_increments; ++col) {
     tmp.item(d_num_increments, col) = 0.0;
   }
   tmp.item(d_num_increments, d_num_increments) = 1.0;

   // d_Up = tmp*A
   delete d_Up;
   d_Up = tmp.mult(A);

   // d_S = sigma.
   delete d_S;
   d_S = sigma;

   // We now have another increment.
   ++d_num_increments;
}

double
IncrementalSVDFastUpdate::checkUOrthogonality()
{
   double result = 0.0;
   if (d_num_increments > 1) {
      int last_col = d_num_increments-1;
      double tmp = 0.0;
      for (int i = 0; i < d_dim; ++i) {
         tmp += d_U->item(i, 0) * d_U->item(i, last_col);
      }
      if (d_size > 1) {
         MPI_Allreduce(&tmp, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
      else {
         result = tmp;
      }
   }
   return result;
}

double
IncrementalSVDFastUpdate::checkUpOrthogonality()
{
   double result = 0.0;
   if (d_num_increments > 1) {
      int last_col = d_num_increments-1;
      for (int i = 0; i < d_num_increments; ++i) {
         result += d_Up->item(i, 0)*d_Up->item(i, last_col);
      }
   }
   return result;
}

void
IncrementalSVDFastUpdate::reOrthogonalize(
   Matrix* m)
{
   CAROM_ASSERT(m != 0);

   int num_rows = m->numRows();
   int num_cols = m->numColumns();
   for (int work = 1; work < num_cols; ++work) {
      double tmp;
      for (int col = 0; col < work; ++col) {
         double factor = 0.0;
         tmp = 0.0;
         for (int i = 0; i < num_rows; ++i) {
            tmp += m->item(i, col)*m->item(i, work);
         }
         if (d_size > 1) {
            MPI_Allreduce(&tmp,
               &factor,
               1,
               MPI_DOUBLE,
               MPI_SUM,
               MPI_COMM_WORLD);
         }
         else {
            factor = tmp;
         }

         for (int i = 0; i < num_rows; ++i) {
            m->item(i, work) -= factor*m->item(i, col);
         }
      }
      double norm = 0.0;
      tmp = 0.0;
      for (int i = 0; i < num_rows; ++i) {
         tmp += m->item(i, work)*m->item(i, work);
      }
      if (d_size > 1) {
         MPI_Allreduce(&tmp, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
      else {
         norm = tmp;
      }
      norm = sqrt(norm);
      for (int i = 0; i < num_rows; ++i) {
         m->item(i, work) /= norm;
      }
   }
}

}
