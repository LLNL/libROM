#include <cmath>
#include <mpi.h>

// Returns the inner product of the 2 distributed vectors v1 and v2.
double
inner_product(
   const double* v1,
   const double* v2,
   int dim,
   int num_procs)
{
   double ip;
   double local_ip = 0.0;
   for (int i = 0; i < dim; ++i) {
      local_ip += v1[i]*v2[i];
   }
   if (num_procs > 1) {
      MPI_Allreduce(&local_ip, &ip, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
   else {
      ip = local_ip;
   }
   return ip;
}

// Returns the norm of the distributed vector v.
double
norm(
   const double* v,
   int dim,
   int num_procs)
{
   double ip = inner_product(v, v, dim, num_procs);
   double norm = sqrt(ip);
   return norm;
}

// Normalizes the distributed vector v and returns its norm.
double
normalize(
   double* v,
   int dim,
   int num_procs)
{
   double Norm = norm(v, dim, num_procs);
   for (int i = 0; i < dim; ++i) {
      v[i] /= Norm;
   }
   return Norm;
}
