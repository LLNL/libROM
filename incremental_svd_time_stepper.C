#include "incremental_svd_time_stepper.h"
#include "incremental_svd_naive.h"
#include "incremental_svd_fast_update.h"

#include <cmath>

namespace CAROM {

incremental_svd_time_stepper::incremental_svd_time_stepper(
   int dim,
   double epsilon,
   bool skip_redundant,
   int increments_per_time_interval,
   double tolerance,
   double max_time_between_increments,
   bool fast_update) :
   d_tol(tolerance),
   d_max_time_between_increments(max_time_between_increments),
   d_next_increment_time(0.0)
{
   if (fast_update) {
      d_isvd.reset(
         new incremental_svd_fast_update(dim,
                                         epsilon,
                                         skip_redundant,
                                         increments_per_time_interval));
   }
   else {
      d_isvd.reset(
         new incremental_svd_naive(dim,
                                   epsilon,
                                   skip_redundant,
                                   increments_per_time_interval));
   }
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
   // Get some preliminary info.
   int dim = d_isvd->getDim();
   int rank = d_isvd->getRank();
   int size = d_isvd->getSize();
   double eps = d_isvd->getEpsilon();

   // Get the current basis vectors.
   const Matrix* basis = getBasis();

   // Compute a bunch of stuff we need.
   Vector u_vec(u_in, dim, true, rank, size);
   Vector* l = basis->TransposeMult(u_vec);
   double cm = u_vec.dot(u_vec);
   double k = cm - (l->dot(*l));
   if (k <= 0) {
      k = 0;
   }
   else {
      k = sqrt(k);
   }

   // Compute j
   Vector* basisl = basis->Mult(*l);
   Vector* j = u_vec.subtract(basisl);
   delete basisl;
   for (int i = 0; i < dim; ++i) {
      j->item(i) /= k;
   }
   delete l;
   delete j;

   Vector rhs_vec(rhs_in, dim, true, rank, size);
   double rhs_norm = rhs_vec.norm();
   double u_norm = sqrt(cm);

   double eps0 = k/(eps+u_norm);
   double deps = (1.0+eps0)*rhs_norm/(1.0e-10+u_norm);
   double dt = (eps-eps0)/(deps+1.0e-10);

   if (dt > d_max_time_between_increments) {
      dt = d_max_time_between_increments;
   }
   else if (dt < 0) {
      dt = 0.0;
   }

   // Return next increment time.
   d_next_increment_time = time + dt;
   return d_next_increment_time;
}

}
