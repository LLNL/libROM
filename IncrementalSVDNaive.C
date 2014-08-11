#include "IncrementalSVDNaive.h"

#include "mpi.h"

#include <cmath>
#include <limits>
#include <stdio.h>

namespace CAROM {

IncrementalSVDNaive::IncrementalSVDNaive(
   int dim,
   double redundancy_tol,
   bool skip_redundant,
   int increments_per_time_interval) :
   IncrementalSVD(dim,
      redundancy_tol,
      skip_redundant,
      increments_per_time_interval),
   d_U(0)
{
}

IncrementalSVDNaive::~IncrementalSVDNaive()
{
   // Delete data members.
   if (d_U) {
      delete d_U;
   }
}

void
IncrementalSVDNaive::increment(
   const double* u_in,
   double time)
{
   // If this is the first SVD then build it.  Otherwise add this increment to
   // the system.
   if (isNewTimeInterval()) {
      buildInitialSVD(u_in, time);
   }
   else {
      buildIncrementalSVD(u_in);
   }

#ifdef DEBUG_ROMS
   const Matrix* basis = getBasis();
   if (d_rank == 0) {
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
      double* m = new double[d_dim*d_num_increments];
      for (int proc = 1; proc < d_size; ++proc) {
         MPI_Status status;
         MPI_Recv(m, d_dim*d_num_increments, MPI_DOUBLE, proc,
                  COMMUNICATE_U, MPI_COMM_WORLD, &status);
         int idx = 0;
         for (int row = 0; row < d_dim; ++row) {
            for (int col = 0; col < d_num_increments; ++col) {
               printf("%.16e ", m[idx++]);
            }
            printf("\n");
         }
      }
      printf("============================================================\n");
      delete [] m;
   }
   else {
      // Send this processor's part of the basis to process 0.
      MPI_Request request;
      MPI_Isend(const_cast<double*>(&basis->item(0, 0)),
                d_dim*d_num_increments, MPI_DOUBLE, 0, COMMUNICATE_U,
                MPI_COMM_WORLD, &request);
   }
#endif
}

const Matrix*
IncrementalSVDNaive::getBasis()
{
   CAROM_ASSERT(d_U);
   d_basis = new Matrix(*d_U);
   return d_basis;
}

void
IncrementalSVDNaive::buildInitialSVD(
   const double* u,
   double time)
{
   // We have a new time interval.

   // If this is not the first time interval then write the basis vectors for
   // the just completed interval.  Delete d_basis and d_S of the just
   // completed time interval.
   if (d_num_time_intervals > 0) {
      if (d_basis) {
         delete d_basis;
         d_basis = 0;
      }
      delete d_U;
      delete d_S;
   }
   ++d_num_time_intervals;
   d_time_interval_start_times.resize(d_num_time_intervals);
   d_time_interval_start_times[d_num_time_intervals-1] = time;

   // Build d_S for this new time interval.
   d_S = new Matrix(1, 1, false, d_rank, d_size);
   Vector u_vec(u, d_dim, true, d_rank, d_size);
   double norm_u = u_vec.norm();
   d_S->item(0, 0) = norm_u;

   // Build d_U for this new time interval.
   d_U = new Matrix(d_dim, 1, true, d_rank, d_size);
   for (int i = 0; i < d_dim; ++i) {
      d_U->item(i, 0) = u[i]/norm_u;
   }

   // We now have the first increment for the new time interval.
   d_num_increments = 1;
}

void
IncrementalSVDNaive::buildIncrementalSVD(
   const double* u)
{
   // l = basis' * u
   Vector u_vec(u, d_dim, true, d_rank, d_size);
   Vector* l = d_U->TransposeMult(u_vec);

   // basisl = basis * l
   Vector* basisl = d_U->Mult(l);

   // Compute k = sqrt(u.u - 2.0*l.l + basisl.basisl) which is ||u - basisl||.
   double k = u_vec.dot(u_vec) - 2.0*l->dot(l) + basisl->dot(basisl);
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
   if (d_epsilon < default_tol) {
      d_epsilon = default_tol;
   }
   bool is_new_increment;
   if (k < d_epsilon) {
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
      delete sigma;
   }
   else if (is_new_increment) {
      // This increment is new.

      // Compute j
      Vector* j = u_vec.subtract(basisl);
      for (int i = 0; i < d_dim; ++i) {
         j->item(i) /= k;
      }

      // addNewIncrement will assign sigma to d_S hence it should not be
      // deleted upon return.
      addNewIncrement(j, A, sigma);
      delete j;

      // Reorthogonalize if necessary.
      int max_U_dim;
      if (d_dim*d_size > d_num_increments) {
         max_U_dim = d_dim*d_size;
      }
      else {
         max_U_dim = d_num_increments;
      }
      if (checkOrthogonality() >
          std::numeric_limits<double>::epsilon()*max_U_dim) {
         reOrthogonalize();
      }
   }
   delete basisl;

   // Clean up.
   delete l;
   delete A;
}

void
IncrementalSVDNaive::addRedundantIncrement(
   const Matrix* A,
   const Matrix* sigma)
{
   // Chop a row and a column off of A to form Amod.  Also form
   // d_S by chopping a row and a column off of sigma.
   Matrix Amod(d_num_increments, d_num_increments, false, d_rank, d_size);
   for (int row = 0; row < d_num_increments; ++row){
      for (int col = 0; col < d_num_increments; ++col) {
         Amod.item(row, col) = A->item(row, col);
         d_S->item(row, col) = sigma->item(row, col);
      }
   }

   // Multiply d_U and Amod and put result into d_U.
   Matrix* U_times_Amod = d_U->Mult(Amod);
   delete d_U;
   d_U = U_times_Amod;
}

void
IncrementalSVDNaive::addNewIncrement(
   const Vector* j,
   const Matrix* A,
   Matrix* sigma)
{
   // Add j as a new column of d_U.  Then multiply by A to form a new d_U.
   Matrix tmp(d_dim, d_num_increments+1, true, d_rank, d_size);
   for (int row = 0; row < d_dim; ++row) {
      for (int col = 0; col < d_num_increments; ++col) {
         tmp.item(row, col) = d_U->item(row, col);
      }
      tmp.item(row, d_num_increments) = j->item(row);
   }
   delete d_U;
   d_U = tmp.Mult(A);

   // d_S = sigma.
   delete d_S;
   d_S = sigma;

   // We now have another increment.
   ++d_num_increments;
}

double
IncrementalSVDNaive::checkOrthogonality()
{
   double result = 0.0;
   if (d_num_increments > 1) {
      int last_col = d_num_increments-1;
      double tmp = 0;
      for (int i = 0; i < d_dim; ++i) {
         tmp += d_U->item(i, 0)*d_U->item(i, last_col);
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

void
IncrementalSVDNaive::reOrthogonalize()
{
   for (int work = 1; work < d_num_increments; ++work) {
      double tmp;
      for (int col = 0; col < work; ++col) {
         double factor = 0.0;
         tmp = 0.0;
         for (int i = 0; i < d_dim; ++i) {
            tmp += d_U->item(i, col)*d_U->item(i, work);
         }
         if (d_size > 1) {
            MPI_Allreduce(&tmp, &factor, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         }
         else {
            factor = tmp;
         }

         for (int i = 0; i < d_dim; ++i) {
            d_U->item(i, work) -= factor*d_U->item(i, col);
         }
      }
      double norm = 0.0;
      tmp = 0.0;
      for (int i = 0; i < d_dim; ++i) {
         tmp += d_U->item(i, work)*d_U->item(i, work);
      }
      if (d_size > 1) {
         MPI_Allreduce(&tmp, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
      else {
         norm = tmp;
      }
      norm = sqrt(norm);
      for (int i = 0; i < d_dim; ++i) {
         d_U->item(i, work) /= norm;
      }
   }
}

}
