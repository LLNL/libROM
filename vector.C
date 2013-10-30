#include "vector.h"
#include <cmath>
#include <mpi.h>

namespace CAROM {

Vector::Vector(
   int dim,
   bool distributed,
   int rank,
   int num_procs) :
   d_dim(dim),
   d_distributed(distributed),
   d_rank(rank),
   d_num_procs(num_procs),
   d_manages_storage(true)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(dim > 0);
   assert(rank < num_procs);
#endif
   d_vec = new double [dim];
}

Vector::Vector(
   double* vec,
   int dim,
   bool distributed,
   int rank,
   int num_procs) :
   d_dim(dim),
   d_distributed(distributed),
   d_rank(rank),
   d_num_procs(num_procs),
   d_manages_storage(false)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(dim > 0);
   assert(rank < num_procs);
#endif
   d_vec = vec;
}

Vector::~Vector()
{
   if (d_manages_storage) {
      delete [] d_vec;
   }
}

double
Vector::dot(const Vector& other) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(d_dim == other.d_dim);
#endif
   double ip;
   double local_ip = 0.0;
   for (int i = 0; i < d_dim; ++i) {
      local_ip += d_vec[i]*other.d_vec[i];
   }
   if (d_num_procs > 1) {
      MPI_Allreduce(&local_ip, &ip, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
   else {
      ip = local_ip;
   }
   return ip;
}

double
Vector::norm() const
{
   double ip = dot(*this);
   double norm = sqrt(ip);
   return norm;
}

double
Vector::normalize()
{
   double Norm = norm();
   for (int i = 0; i < d_dim; ++i) {
      d_vec[i] /= Norm;
   }
   return Norm;
}

}
