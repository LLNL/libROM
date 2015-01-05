/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The class that determines the next time at which a sample
 *              should be taken for basis generation using an incremental SVD
 *              approach.
 *
 *****************************************************************************/

#include "IncrementalSVDSampler.h"
#include "IncrementalSVDNaive.h"
#include "IncrementalSVDFastUpdate.h"

#include <cmath>

namespace CAROM {

IncrementalSVDSampler::IncrementalSVDSampler(
   int dim,
   double redundancy_tol,
   bool skip_redundant,
   int samples_per_time_interval,
   double sampling_tol,
   double max_time_between_samples,
   bool fast_update,
   bool debug_rom) :
   d_tol(sampling_tol),
   d_max_time_between_samples(max_time_between_samples),
   d_next_sample_time(0.0)
{
   CAROM_ASSERT(dim > 0);
   CAROM_ASSERT(redundancy_tol > 0.0);
   CAROM_ASSERT(samples_per_time_interval > 0);
   CAROM_ASSERT(sampling_tol > 0.0);
   CAROM_ASSERT(max_time_between_samples > 0.0);

   if (fast_update) {
      d_isvd.reset(
         new IncrementalSVDFastUpdate(dim,
            redundancy_tol,
            skip_redundant,
            samples_per_time_interval,
            debug_rom));
   }
   else {
      d_isvd.reset(
         new IncrementalSVDNaive(dim,
            redundancy_tol,
            skip_redundant,
            samples_per_time_interval,
            debug_rom));
   }
}

IncrementalSVDSampler::~IncrementalSVDSampler()
{
}

double
IncrementalSVDSampler::computeNextSampleTime(
   double* u_in,
   double* rhs_in,
   double time)
{
   CAROM_ASSERT(u_in != 0);
   CAROM_ASSERT(rhs_in != 0);
   CAROM_ASSERT(time >= 0.0);

   // Get some preliminary info.
   int dim = d_isvd->getDim();

   // Get the current basis vectors.
   const Matrix* basis = getBasis();

   // Compute the square of the norm of the projection of the normalized
   // current state into the reduced order space.
   Vector u_vec(u_in, dim, true);
   u_vec.normalize();
   Vector* l = basis->transposeMult(u_vec);

   // Compute the square of the norm of the normalized current state.
   double cm = u_vec.inner_product(u_vec);

   // Compute the norm of the out of plane component.
   double proj_error_norm = cm - (l->inner_product(l));
   if (proj_error_norm <= 0) {
      proj_error_norm = 0.0;
   }
   else {
      proj_error_norm = sqrt(proj_error_norm);
   }
   delete l;

   // Compute the square of the norm of the projection of the normalized rhs
   // into the reduced order space.
   Vector rhs_vec(rhs_in, dim, true);
   rhs_vec.normalize();
   l = basis->transposeMult(rhs_vec);

   // Compute the square of the norm of the normalized rhs.
   cm = rhs_vec.inner_product(rhs_vec);

   // Compute the norm of the out of plane component.
   double proj_error_deriv_norm = cm - (l->inner_product(l));
   if (proj_error_deriv_norm <= 0) {
      proj_error_deriv_norm = 0.0;
   }
   else {
      proj_error_deriv_norm = sqrt(proj_error_deriv_norm);
   }
   delete l;

   // Compute dt from these two norms.
   double dt;
   if (proj_error_deriv_norm > 0) {
      dt = (d_tol - proj_error_norm) / proj_error_deriv_norm;
      if (dt > d_max_time_between_samples) {
         dt = d_max_time_between_samples;
      }
      else if (dt < 0) {
         dt = 0.0;
      }
   }
   else {
      dt = d_max_time_between_samples;
   }

   // Return next sample time.
   d_next_sample_time = time + dt;
   return d_next_sample_time;
}

}
