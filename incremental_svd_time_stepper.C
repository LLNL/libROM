#include "incremental_svd_time_stepper.h"
#include <cmath>

namespace CAROM {

incremental_svd_time_stepper::incremental_svd_time_stepper(
   int* argc,
   char*** argv,
   int dim,
   double epsilon,
   bool skip_redundant,
   int max_time_steps_between_increments) :
   d_max_time_steps_between_increments(max_time_steps_between_increments),
   d_isvd(new incremental_svd(argc,
                              argv,
                              dim,
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
   int rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   int dim = d_isvd->getDim();

   // Convert representation free C arrays of doubles into a PETSc vectors.
   Vec u;
   VecCreate(PETSC_COMM_WORLD, &u);
   VecSetSizes(u, dim, PETSC_DECIDE);
   VecSetType(u, VECMPI);
   int offset = rank*dim;
   int* vec_locs = new int [dim];
   for (int i = 0; i < dim; ++i) {
      vec_locs[i] = offset + i;
   }
   VecSetValues(u, dim, vec_locs, u_in, INSERT_VALUES);
   VecAssemblyBegin(u);
   VecAssemblyEnd(u);
   VecView(u, PETSC_VIEWER_STDOUT_WORLD);

   Vec rhs;
   VecCreate(PETSC_COMM_WORLD, &rhs);
   VecSetSizes(rhs, dim, PETSC_DECIDE);
   VecSetType(rhs, VECMPI);
   VecSetValues(rhs, dim, vec_locs, rhs_in, INSERT_VALUES);
   VecAssemblyBegin(rhs);
   VecAssemblyEnd(rhs);
   VecView(rhs, PETSC_VIEWER_STDOUT_WORLD);
   delete [] vec_locs;

   // Get the norm of J from the incremental svd algorithm.
   double norm_j = d_isvd->computeNormJ(u);

   // Compute norm_u = sqrt(u.u).
   double norm_u;
   VecDot(u, u, &norm_u);
   norm_u = sqrt(norm_u);

   // Compute the norm_rhs = sqrt(rhs.rhs).
   double norm_rhs;
   VecDot(rhs, rhs, &norm_rhs);
   norm_rhs = sqrt(norm_rhs);

   // Compute delta t to next increment time.
   double epsilon = d_isvd->getEpsilon();
   double eps0 = norm_j/(epsilon+norm_u);
   double deps = (1.0 + eps0)*norm_rhs/(1.0e-10 + norm_u);
   double dtcheck = (epsilon - eps0)/(deps + 1e-10);

   // Clean up.    
   VecDestroy(&u);
   VecDestroy(&rhs);

   // Return next increment time.
   return time + dtcheck;
}

}
