/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The concrete implementation of the incremental SVD algorithm
 *              that is equivalent to but computationally more expensive than
 *              the "fast update" method.
 *
 *****************************************************************************/

#include "IncrementalSVDNaive.h"

#include "mpi.h"

#include <cmath>
#include <limits>

namespace CAROM {

IncrementalSVDNaive::IncrementalSVDNaive(
   int dim,
   double redundancy_tol,
   bool skip_redundant,
   int samples_per_time_interval,
   bool debug_rom) :
   IncrementalSVD(dim,
      redundancy_tol,
      skip_redundant,
      samples_per_time_interval,
      debug_rom),
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
   CAROM_ASSERT(u != 0);
   CAROM_ASSERT(time >= 0.0);

   // We have a new time interval.

   // If this is not the first time interval then write the basis vectors for
   // the just completed interval.  Delete d_basis and d_S of the just
   // completed time interval.
   int num_time_intervals =
      static_cast<int>(d_time_interval_start_times.size());
   if (num_time_intervals > 0) {
      if (d_basis) {
         delete d_basis;
         d_basis = 0;
      }
      delete d_U;
      delete d_S;
   }
   d_time_interval_start_times.resize(num_time_intervals+1);
   d_time_interval_start_times[num_time_intervals] = time;

   // Build d_S for this new time interval.
   d_S = new Matrix(1, 1, false);
   Vector u_vec(u, d_dim, true);
   double norm_u = u_vec.norm();
   d_S->item(0, 0) = norm_u;

   // Build d_U for this new time interval.
   d_U = new Matrix(d_dim, 1, true);
   for (int i = 0; i < d_dim; ++i) {
      d_U->item(i, 0) = u[i]/norm_u;
   }

   // We now have the first sample for the new time interval.
   d_num_samples = 1;
}

void
IncrementalSVDNaive::buildIncrementalSVD(
   const double* u)
{
   CAROM_ASSERT(u != 0);

   // l = basis' * u
   Vector u_vec(u, d_dim, true);
   Vector* l = d_U->transposeMult(u_vec);

   // basisl = basis * l
   Vector* basisl = d_U->mult(l);

   // Compute k = sqrt(u.u - 2.0*l.l + basisl.basisl) which is ||u - basisl||.
   double k = u_vec.inner_product(u_vec) - 2.0*l->inner_product(l) +
      basisl->inner_product(basisl);
   if (k <= 0) {
      k = 0;
   }
   else {
      k = sqrt(k);
   }

   // Use k to see if this sample is new.
   double default_tol;
   if (d_num_samples == 1) {
      default_tol = std::numeric_limits<double>::epsilon();
   }
   else {
      default_tol = d_num_samples*std::numeric_limits<double>::epsilon()*d_S->item(0, 0);
   }
   if (d_redundancy_tol < default_tol) {
      d_redundancy_tol = default_tol;
   }
   bool is_new_sample;
   if (k < d_redundancy_tol) {
      k = 0;
      is_new_sample = false;
   }
   else {
      is_new_sample = true;
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

   // We need to add the sample if it is new or if it is redundant and we are
   // not skipping redundant samples.
   if (!is_new_sample && !d_skip_redundant) {
      // This sample is redundant and we are not skipping redundant samples.
      addRedundantSample(A, sigma);
      delete sigma;
   }
   else if (is_new_sample) {
      // This sample is new.

      // Compute j
      Vector* j = u_vec.minus(basisl);
      for (int i = 0; i < d_dim; ++i) {
         j->item(i) /= k;
      }

      // addNewSample will assign sigma to d_S hence it should not be deleted
      // upon return.
      addNewSample(j, A, sigma);
      delete j;

      // Reorthogonalize if necessary.
      long int max_U_dim;
      if (d_total_dim > d_num_samples) {
         max_U_dim = d_total_dim;
      }
      else {
         max_U_dim = d_num_samples;
      }
      if (checkOrthogonality() >
          std::numeric_limits<double>::epsilon()*max_U_dim) {
         reOrthogonalize(d_U);
      }
   }
   delete basisl;

   // Clean up.
   delete l;
   delete A;
}

void
IncrementalSVDNaive::addRedundantSample(
   const Matrix* A,
   const Matrix* sigma)
{
   CAROM_ASSERT(A != 0);
   CAROM_ASSERT(sigma != 0);

   // Chop a row and a column off of A to form Amod.  Also form
   // d_S by chopping a row and a column off of sigma.
   Matrix Amod(d_num_samples, d_num_samples, false);
   for (int row = 0; row < d_num_samples; ++row){
      for (int col = 0; col < d_num_samples; ++col) {
         Amod.item(row, col) = A->item(row, col);
         d_S->item(row, col) = sigma->item(row, col);
      }
   }

   // Multiply d_U and Amod and put result into d_U.
   Matrix* U_times_Amod = d_U->mult(Amod);
   delete d_U;
   d_U = U_times_Amod;
}

void
IncrementalSVDNaive::addNewSample(
   const Vector* j,
   const Matrix* A,
   Matrix* sigma)
{
   // Add j as a new column of d_U.  Then multiply by A to form a new d_U.
   Matrix tmp(d_dim, d_num_samples+1, true);
   for (int row = 0; row < d_dim; ++row) {
      for (int col = 0; col < d_num_samples; ++col) {
         tmp.item(row, col) = d_U->item(row, col);
      }
      tmp.item(row, d_num_samples) = j->item(row);
   }
   delete d_U;
   d_U = tmp.mult(A);

   // d_S = sigma.
   delete d_S;
   d_S = sigma;

   // We now have another sample.
   ++d_num_samples;
}

double
IncrementalSVDNaive::checkOrthogonality()
{
   double result = 0.0;
   if (d_num_samples > 1) {
      int last_col = d_num_samples-1;
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

}
