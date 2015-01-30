/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2015 Lawrence Livermore National Security, LLC
 * Description: A simple, parallel Vector class with the utility needed to
 *              support the basis generation methods of this library.  A
 *              distributed Vector has its rows distributed across processors.
 *
 *****************************************************************************/

#include "Vector.h"

#include "mpi.h"

#include <cmath>
#include <string.h>

namespace CAROM {

Vector::Vector(
   int dim,
   bool distributed) :
   d_dim(dim),
   d_distributed(distributed)
{
   CAROM_ASSERT(dim > 0);
   d_vec = new double [dim];
}

Vector::Vector(
   const double* vec,
   int dim,
   bool distributed) :
   d_dim(dim),
   d_distributed(distributed)
{
   CAROM_ASSERT(vec != 0);
   CAROM_ASSERT(dim > 0);
   d_vec = new double [dim];
   memcpy(d_vec, vec, dim*sizeof(double));
}

Vector::Vector(
   const Vector& other) :
   d_dim(other.d_dim),
   d_distributed(other.d_distributed)
{
   d_vec = new double [d_dim];
   memcpy(d_vec, other.d_vec, d_dim*sizeof(double));
}

Vector::~Vector()
{
   delete [] d_vec;
}

Vector&
Vector::operator = (
   const Vector& rhs)
{
   d_distributed = rhs.d_distributed;
   if (d_dim != rhs.d_dim) {
      d_dim = rhs.d_dim;
      delete [] d_vec;
      d_vec = new double[d_dim];
   }
   memcpy(d_vec, rhs.d_vec, d_dim*sizeof(double));
   return *this;
}

double
Vector::inner_product(
   const Vector& other) const
{
   CAROM_ASSERT(dim() == other.dim());
   CAROM_ASSERT(distributed() == other.distributed());
   double ip;
   double local_ip = 0.0;
   for (int i = 0; i < d_dim; ++i) {
      local_ip += d_vec[i]*other.d_vec[i];
   }
   int mpi_init;
   MPI_Initialized(&mpi_init);
   int num_procs;
   if (mpi_init) {
      MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   }
   else {
      num_procs = 1;
   }
   if (num_procs > 1 && d_distributed) {
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
   double ip = inner_product(this);
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
Vector::plus(
   const Vector& other) const
{
   CAROM_ASSERT(distributed() == other.distributed());
   CAROM_ASSERT(dim() == other.dim());

   Vector* result = new Vector(d_dim, d_distributed);
   for (int i = 0; i < d_dim; ++i) {
      result->d_vec[i] = d_vec[i] + other.d_vec[i];
   }
   return result;
}

Vector*
Vector::minus(
   const Vector& other) const
{
   CAROM_ASSERT(distributed() == other.distributed());
   CAROM_ASSERT(dim() == other.dim());

   Vector* result = new Vector(d_dim, d_distributed);
   for (int i = 0; i < d_dim; ++i) {
      result->d_vec[i] = d_vec[i] - other.d_vec[i];
   }
   return result;
}

}
