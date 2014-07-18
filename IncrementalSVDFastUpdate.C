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
   int increments_per_time_interval) :
   IncrementalSVD(dim,
      redundancy_tol,
      skip_redundant,
      increments_per_time_interval),
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
   // We have a new time interval.

   // If this is not the first time interval then delete d_basis, d_U, d_Up,
   // and d_S of the just completed time interval.
   if (d_num_time_intervals > 0) {
      delete d_basis;
      delete d_U;
      delete d_Up;
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

   // Build d_Up for this new time interval.
   d_Up = new Matrix(1, 1, false, d_rank, d_size);
   d_Up->item(0, 0) = 1.0;

   // Build d_U for this new time interval.
   d_U = new Matrix(d_dim, 1, true, d_rank, d_size);
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
   // l = basis' * u
   Vector u_vec(u, d_dim, true, d_rank, d_size);
   Vector* l = d_basis->TransposeMult(u_vec);

   // basisl = basis * l
   Vector* basisl = d_basis->Mult(*l);

   // Compute k = sqrt(u.u - 2.0*l.l + basisl.basisl) which is ||u - basisl||.
   double k = u_vec.dot(u_vec) - 2.0*l->dot(*l) + basisl->dot(*basisl);
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

      // Reorthogonalize if necessary.
      if (fabs(checkUpOrthogonality()) >
          std::numeric_limits<double>::epsilon()*d_num_increments) {
         reOrthogonalizeUp();
      }

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
      if (fabs(checkUOrthogonality()) >
          std::numeric_limits<double>::epsilon()*max_U_dim) {
         reOrthogonalizeU();
      }
      if (fabs(checkUpOrthogonality()) >
          std::numeric_limits<double>::epsilon()*d_num_increments) {
         reOrthogonalizeUp();
      }
   }
   delete basisl;

   // Compute the basis vectors.
   delete d_basis;
   d_basis = d_U->Mult(*d_Up);

   // Clean up.
   delete l;
   delete A;
}

void
IncrementalSVDFastUpdate::addRedundantIncrement(
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

   // Multiply d_Up and Amod and put result into d_Up.
   Matrix* Up_times_Amod = d_Up->Mult(Amod);
   delete d_Up;
   d_Up = Up_times_Amod;
}

void
IncrementalSVDFastUpdate::addNewIncrement(
   const Vector* j,
   const Matrix* A,
   Matrix* sigma)
{
   // Add j as a new column of d_U.
   Matrix* newU = new Matrix(d_dim, d_num_increments+1, true, d_rank, d_size);
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
   Matrix tmp(d_num_increments+1, d_num_increments+1, false, d_rank, d_size);
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
   d_Up = tmp.Mult(*A);

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
IncrementalSVDFastUpdate::reOrthogonalizeU()
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

void
IncrementalSVDFastUpdate::reOrthogonalizeUp()
{
   for (int work = 1; work < d_num_increments; ++work) {
      double tmp;
      for (int col = 0; col < work; ++col) {
         double factor = 0.0;
         tmp = 0.0;
         for (int i = 0; i < d_num_increments; ++i) {
            tmp += d_Up->item(i, col)*d_Up->item(i, work);
         }
         if (d_size > 1) {
            MPI_Allreduce(&tmp, &factor, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         }
         else {
            factor = tmp;
         }

         for (int i = 0; i < d_num_increments; ++i) {
            d_Up->item(i, work) -= factor*d_Up->item(i, col);
         }
      }
      double norm = 0.0;
      tmp = 0.0;
      for (int i = 0; i < d_num_increments; ++i) {
         tmp += d_Up->item(i, work)*d_Up->item(i, work);
      }
      if (d_size > 1) {
         MPI_Allreduce(&tmp, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
      else {
         norm = tmp;
      }
      norm = sqrt(norm);
      for (int i = 0; i < d_num_increments; ++i) {
         d_Up->item(i, work) /= norm;
      }
   }
}

}
