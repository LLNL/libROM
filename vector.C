#include "vector.h"

#include "mpi.h"

#include <cmath>
#include <string.h>

namespace CAROM {

Vector::Vector(
   int dim,
   bool distributed,
   int rank,
   int num_procs) :
   d_dim(dim),
   d_distributed(distributed),
   d_rank(rank),
   d_num_procs(num_procs)
{
   CAROM_ASSERT(dim > 0);
   CAROM_ASSERT(rank < num_procs);
   d_vec = new double [dim];
}

Vector::Vector(
   const double* vec,
   int dim,
   bool distributed,
   int rank,
   int num_procs) :
   d_dim(dim),
   d_distributed(distributed),
   d_rank(rank),
   d_num_procs(num_procs)
{
   CAROM_ASSERT(dim > 0);
   CAROM_ASSERT(rank < num_procs);
   d_vec = new double [dim];
   memcpy(d_vec, vec, dim*sizeof(double));
}

Vector::~Vector()
{
   delete [] d_vec;
}

double
Vector::dot(
   const Vector& other) const
{
   CAROM_ASSERT(d_dim == other.d_dim);
   double ip;
   double local_ip = 0.0;
   for (int i = 0; i < d_dim; ++i) {
      local_ip += d_vec[i]*other.d_vec[i];
   }
   if (d_num_procs > 1 && d_distributed) {
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

Vector*
Vector::add(const Vector* other) const
{
   CAROM_ASSERT(d_distributed == other->d_distributed);
   CAROM_ASSERT(d_dim == other->d_dim);

   Vector* result = new Vector(d_dim, d_distributed, d_rank, d_num_procs);
   for (int i = 0; i < d_dim; ++i) {
      result->d_vec[i] = d_vec[i] + other->d_vec[i];
   }
   return result;
}

Vector*
Vector::subtract(const Vector* other) const
{
   CAROM_ASSERT(d_distributed == other->d_distributed);
   CAROM_ASSERT(d_dim == other->d_dim);

   Vector* result = new Vector(d_dim, d_distributed, d_rank, d_num_procs);
   for (int i = 0; i < d_dim; ++i) {
      result->d_vec[i] = d_vec[i] - other->d_vec[i];
   }
   return result;
}

}
