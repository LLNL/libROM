#include "incremental_svd_time_stepper.h"
#include "vector_utils.h"
#include <cmath>
#include <mpi.h>

namespace CAROM {

incremental_svd_time_stepper::incremental_svd_time_stepper(
   int dim,
   double epsilon,
   bool skip_redundant,
   int increments_per_time_interval,
   int max_time_steps_between_increments) :
   d_max_time_steps_between_increments(max_time_steps_between_increments),
   d_next_increment_time(0.0),
   d_isvd(new incremental_svd(dim,
                              epsilon,
                              skip_redundant,
                              increments_per_time_interval))
{
}

incremental_svd_time_stepper::~incremental_svd_time_stepper()
{
}

double
incremental_svd_time_stepper::computeNextIncrementTime(
   double* u_in,
   double* rhs_in,
   double time)
{
   int dim = d_isvd->getDim();
   int num_procs = d_isvd->getSize();

   // Get the norm of J from the incremental svd algorithm.
   double norm_j = d_isvd->getNormJ();

   // Compute the norm of u_in.
   double norm_u = norm(u_in, dim, num_procs);

   // Compute the norm of rhs_in.
   double norm_rhs = norm(rhs_in, dim, num_procs);

   // Compute delta t to next increment time.
   double epsilon = d_isvd->getEpsilon();
   double eps0 = norm_j/(epsilon+norm_u);
   double deps = (1.0 + eps0)*norm_rhs/(1.0e-10 + norm_u);
   double dtcheck = (epsilon - eps0)/(deps + 1e-10);

   // Return next increment time.
   d_next_increment_time = time + dtcheck;
   return d_next_increment_time;
}

}
