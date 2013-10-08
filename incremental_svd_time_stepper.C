#include "incremental_svd_time_stepper.h"
#include <cmath>
#include <mpi.h>

namespace CAROM {

incremental_svd_time_stepper::incremental_svd_time_stepper(
   int dim,
   double epsilon,
   bool skip_redundant,
   int max_time_steps_between_increments) :
   d_max_time_steps_between_increments(max_time_steps_between_increments),
   d_next_increment_time(0.0),
   d_isvd(new incremental_svd(dim,
                              epsilon,
                              skip_redundant))
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
   int mpi_init;
   MPI_Initialized(&mpi_init);
   int dim = d_isvd->getDim();

   // Get the norm of J from the incremental svd algorithm.
   double norm_j = d_isvd->getNormJ();

   // Compute the norm of u.
   double norm_u = 0;
   double tmp = 0;
   for (int i = 0; i < dim; ++i) {
      tmp += u_in[i]*u_in[i];
   }

   if (mpi_init) {
      MPI_Allreduce(&tmp, &norm_u, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
   else {
      norm_u = tmp;
   }
   norm_u = sqrt(norm_u);

   // Compute the norm of rhs.
   double norm_rhs = 0.0;
   tmp = 0.0;
   for (int i = 0; i < dim; ++i) {
      tmp += rhs_in[i]*rhs_in[i];
   }
   if (mpi_init) {
      MPI_Allreduce(&tmp, &norm_rhs, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
   else {
      norm_rhs = tmp;
   }
   norm_rhs = sqrt(norm_rhs);

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
