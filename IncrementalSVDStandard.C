/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The concrete implementation of the incremental SVD algorithm
//              that is equivalent to but computationally more expensive than
//              the "fast update" method.

#include "IncrementalSVDStandard.h"
#include "HDFDatabase.h"

#include "mpi.h"

#include <cmath>
#include <limits>

namespace CAROM {

IncrementalSVDStandard::IncrementalSVDStandard(
   int dim,
   double linearity_tol,
   bool skip_linearly_dependent,
   int max_basis_dimension,
   int samples_per_time_interval,
   const std::string& basis_file_name,
   bool save_state,
   bool restore_state,
   bool updateRightSV,
   bool debug_algorithm) :
   IncrementalSVD(dim,
      linearity_tol,
      skip_linearly_dependent,
      max_basis_dimension,
      samples_per_time_interval,
      basis_file_name,
      save_state,
      restore_state,
      updateRightSV,
      debug_algorithm)
{
   CAROM_ASSERT(dim > 0);
   CAROM_ASSERT(linearity_tol > 0.0);
   CAROM_ASSERT(samples_per_time_interval > 0);

   // If the state of the SVD is to be restored, do it now.  The base class,
   // IncrementalSVD, has already opened the database and restored the state
   // common to all incremental algorithms.  This particular class has no other
   // state to read and only needs to compute the basis.  If the database could
   // not be found then we can not restore the state.
   if (restore_state && d_state_database) {
      // Close and delete the database.
      d_state_database->close();
      delete d_state_database;

      // Compute the basis.
      computeBasis();
   }
}

IncrementalSVDStandard::~IncrementalSVDStandard()
{
   // If the state of the SVD is to be saved, then create the database now.
   // The IncrementalSVD base class destructor will save d_S and d_U.  This
   // derived class has no specific state to save.
   //
   // If there are multiple time intervals then saving and restoring the state
   // does not make sense as there is not one, all encompassing, basis.
   if (d_save_state && d_time_interval_start_times.size() == 1) {
      // Create state database file.
      d_state_database = new HDFDatabase();
      d_state_database->create(d_state_file_name);
   }
}

void
IncrementalSVDStandard::buildInitialSVD(
   double* u,
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
      delete d_basis;
      delete d_U;
      delete d_S;
      delete d_W;
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

   // Build d_W for this new time interval.
   if (d_updateRightSV) {
     d_W = new Matrix(1, 1, false);
     d_W->item(0, 0) = 1.0;
   }

   // Compute the basis vectors for this time interval.
   computeBasis();

   // We now have the first sample for the new time interval.
   d_num_samples = 1;
}

void
IncrementalSVDStandard::computeBasis()
{
   /* Invalidate existing cached basis and update cached basis */
   delete d_basis;
   d_basis = new Matrix(*d_U);

   if (d_updateRightSV)
   {
      delete d_basis_right;
      d_basis_right = new Matrix(*d_W);
   }
}

void
IncrementalSVDStandard::addLinearlyDependentSample(
   const Matrix* A,
   const Matrix* W,
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

   // Chop a column off of W to form Wmod.
   Matrix* new_d_W;
   if (d_updateRightSV) {
     new_d_W = new Matrix(d_num_rows_of_W+1, d_num_samples, false);
     for (int row = 0; row < d_num_rows_of_W; ++row) {
        for (int col = 0; col < d_num_samples; ++col) {
           double new_d_W_entry = 0.0;
           for (int entry = 0; entry < d_num_samples; ++entry) {
              new_d_W_entry += d_W->item(row, entry)*W->item(entry, col);
           }
           new_d_W->item(row, col) = new_d_W_entry;
        }
     }
     for (int col = 0; col < d_num_samples; ++col) {
        new_d_W->item(d_num_rows_of_W, col) = W->item(d_num_samples, col);
     }
     delete d_W;
     d_W = new_d_W;
     ++d_num_rows_of_W;
   }

   // Reorthogonalize if necessary.
   long int max_U_dim;
   if (d_num_samples > d_total_dim) {
      max_U_dim = d_num_samples;
   }
   else {
      max_U_dim = d_total_dim;
   }
   if (fabs(checkOrthogonality(d_U)) >
       std::numeric_limits<double>::epsilon()*static_cast<double>(max_U_dim)) {
      reOrthogonalize(d_U);
   }
}

void
IncrementalSVDStandard::addNewSample(
   const Vector* j,
   const Matrix* A,
   const Matrix* W,
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

   Matrix* new_d_W;
   if (d_updateRightSV) {
     new_d_W = new Matrix(d_num_rows_of_W+1, d_num_samples+1, false);
     for (int row = 0; row < d_num_rows_of_W; ++row) {
        for (int col = 0; col < d_num_samples+1; ++col) {
           double new_d_W_entry = 0.0;
           for (int entry = 0; entry < d_num_samples; ++entry) {
              new_d_W_entry += d_W->item(row, entry)*W->item(entry, col);
           }
           new_d_W->item(row, col) = new_d_W_entry;
        }
     }
     for (int col = 0; col < d_num_samples+1; ++col) {
        new_d_W->item(d_num_rows_of_W, col) = W->item(d_num_samples, col);
     }
     delete d_W;
     d_W = new_d_W;
   }

   delete d_S;
   d_S = sigma;

   // We now have another sample.
   ++d_num_samples;
   ++d_num_rows_of_W;

   // Reorthogonalize if necessary.
   long int max_U_dim;
   if (d_num_samples > d_total_dim) {
      max_U_dim = d_num_samples;
   }
   else {
      max_U_dim = d_total_dim;
   }
   if (fabs(checkOrthogonality(d_U)) >
       std::numeric_limits<double>::epsilon()*static_cast<double>(max_U_dim)) {
      reOrthogonalize(d_U);
   }
   if (d_updateRightSV) {
     if (fabs(checkOrthogonality(d_W)) >
         std::numeric_limits<double>::epsilon()*d_num_samples) {
        reOrthogonalize(d_W);
     }
   }
}

}
