#include "IncrementalSVDTimeStepper.h"
#include "IncrementalSVDNaive.h"
#include "IncrementalSVDFastUpdate.h"

#include <cmath>

namespace CAROM {

IncrementalSVDTimeStepper::IncrementalSVDTimeStepper(
   int dim,
   double redundancy_tol,
   bool skip_redundant,
   int increments_per_time_interval,
   double sampling_tol,
   double max_time_between_increments,
   bool fast_update) :
   d_tol(sampling_tol),
   d_max_time_between_increments(max_time_between_increments),
   d_next_increment_time(0.0)
{
   CAROM_ASSERT(sampling_tol > 0.0);

   if (fast_update) {
      d_isvd.reset(
         new IncrementalSVDFastUpdate(dim,
                                      redundancy_tol,
                                      skip_redundant,
                                      increments_per_time_interval));
   }
   else {
      d_isvd.reset(
         new IncrementalSVDNaive(dim,
                                 redundancy_tol,
                                 skip_redundant,
                                 increments_per_time_interval));
   }
}

IncrementalSVDTimeStepper::~IncrementalSVDTimeStepper()
{
}

double
IncrementalSVDTimeStepper::computeNextIncrementTime(
   double* u_in,
   double* rhs_in,
   double time)
{
   // Get some preliminary info.
   int dim = d_isvd->getDim();
   int rank = d_isvd->getRank();
   int size = d_isvd->getSize();

   // Get the current basis vectors.
   const Matrix* basis = getBasis();

   // Compute the projection of the current state into the reduced order space.
   Vector u_vec(u_in, dim, true, rank, size);
   Vector* l = basis->TransposeMult(u_vec);

   // Compute the difference in the norms of the current state and its
   // projection into the reduced order space.
   double cm = u_vec.dot(u_vec);
   double proj_error_norm = cm - (l->dot(*l));
   if (proj_error_norm <= 0) {
      proj_error_norm = 0.0;
   }
   else {
      proj_error_norm = sqrt(proj_error_norm);
   }
   delete l;

   // Compute the projection of the rhs into the reduced order space.
   Vector rhs_vec(rhs_in, dim, true, rank, size);
   l = basis->TransposeMult(rhs_vec);

   // Compute the difference in the norms of the rhs and its projection into
   // the reduced order space.
   cm = rhs_vec.dot(rhs_vec);
   double proj_error_deriv_norm = cm - (l->dot(*l));
   if (proj_error_deriv_norm <= 0) {
      proj_error_deriv_norm = 0.0;
   }
   else {
      proj_error_deriv_norm = sqrt(proj_error_deriv_norm);
   }
   delete l;

   // Compute dt from these two norms.
   double dt;
   if (proj_error_norm > 0) {
      dt = (d_tol - proj_error_norm) / proj_error_deriv_norm;
      if (dt > d_max_time_between_increments) {
         dt = d_max_time_between_increments;
      }
      else if (dt < 0) {
         dt = 0.0;
      }
   }
   else {
      dt = d_max_time_between_increments;
   }

   // Return next increment time.
   d_next_increment_time = time + dt;
   return d_next_increment_time;
}

}
