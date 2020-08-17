/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: A simple, parallel Vector class with the utility needed to
//              support the basis generation methods of this library.  A
//              distributed Vector has its rows distributed across processors.

#include "Vector.h"
#include "HDFDatabase.h"

#include "mpi.h"

#include <cmath>
#include <string.h>

namespace CAROM {

Vector::Vector() :
   d_vec(NULL),
   d_alloc_size(0),
   d_distributed(false),
   d_owns_data(true)
{}

Vector::Vector(
   int dim,
   bool distributed) :
   d_vec(NULL),
   d_alloc_size(0),
   d_distributed(distributed),
   d_owns_data(true)
{
   CAROM_ASSERT(dim > 0);
   setSize(dim);
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init) {
      MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
   }
   else {
      d_num_procs = 1;
   }
}

Vector::Vector(
   double* vec,
   int dim,
   bool distributed,
   bool copy_data) :
   d_vec(NULL),
   d_alloc_size(0),
   d_distributed(distributed),
   d_owns_data(copy_data)
{
   CAROM_ASSERT(vec != 0);
   CAROM_ASSERT(dim > 0);
   if (copy_data) {
      setSize(dim);
      memcpy(d_vec, vec, d_alloc_size*sizeof(double));
   }
   else {
      d_vec = vec;
      d_alloc_size = dim;
      d_dim = dim;
   }
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init) {
      MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
   }
   else {
      d_num_procs = 1;
   }
}

Vector::Vector(
   const Vector& other) :
   d_vec(NULL),
   d_alloc_size(0),
   d_distributed(other.d_distributed),
   d_owns_data(true)
{
   setSize(other.d_dim);
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init) {
      MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
   }
   else {
      d_num_procs = 1;
   }
   memcpy(d_vec, other.d_vec, d_alloc_size*sizeof(double));
}

Vector::~Vector()
{
   if (d_owns_data && d_vec) {
      delete [] d_vec;
   }
}

Vector&
Vector::operator = (
   const Vector& rhs)
{
   d_distributed = rhs.d_distributed;
   d_num_procs = rhs.d_num_procs;
   setSize(rhs.d_dim);
   memcpy(d_vec, rhs.d_vec, d_dim*sizeof(double));
   return *this;
}

Vector&
Vector::operator += (
   const Vector& rhs)
{
   CAROM_ASSERT(d_dim == rhs.d_dim);
   for(int i=0; i<d_dim; ++i) d_vec[i] += rhs.d_vec[i];
   return *this;
}

Vector&
Vector::operator = (const double& a)
{
   for(int i=0; i<d_dim; ++i) d_vec[i] = a;
   return *this;
}

Vector&
Vector::transform(void (*f) (const int size, double* vector)) {
      (*f)(d_dim, d_vec);
      return *this;
}

void
Vector::transform(Vector& result, void (*f) (const int size, double* vector)) {
      result = *this;
      (*f)(d_dim, result.getVector());
}

void
Vector::transform(Vector*& result, void (*f) (const int size, double* vector)) {
      // If the result has not been allocated then do so.  Otherwise size it
      // correctly.
      if (result == 0) {
         result = new Vector(d_dim, d_distributed);
      }
      *result = *this;
      (*f)(d_dim, result->getVector());
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

void
Vector::plus(
   const Vector& other,
   Vector*& result) const
{
   CAROM_ASSERT(result == 0 || result->distributed() == distributed());
   CAROM_ASSERT(distributed() == other.distributed());
   CAROM_ASSERT(dim() == other.dim());

   // If the result has not been allocated then do so.  Otherwise size it
   // correctly.
   if (result == 0) {
      result = new Vector(d_dim, d_distributed);
   }
   else {
      result->setSize(d_dim);
   }

   // Do the addition.
   for (int i = 0; i < d_dim; ++i) {
      result->d_vec[i] = d_vec[i] + other.d_vec[i];
   }
}

void
Vector::plus(
   const Vector& other,
   Vector& result) const
{
   CAROM_ASSERT(result.distributed() == distributed());
   CAROM_ASSERT(distributed() == other.distributed());
   CAROM_ASSERT(dim() == other.dim());

   // Size result correctly.
   result.setSize(d_dim);

   // Do the addition.
   for (int i = 0; i < d_dim; ++i) {
      result.d_vec[i] = d_vec[i] + other.d_vec[i];
   }
}

void
Vector::plusAx(
   double factor,
   const Vector& other,
   Vector*& result) const
{
   CAROM_ASSERT(result == 0 || result->distributed() == distributed());
   CAROM_ASSERT(distributed() == other.distributed());
   CAROM_ASSERT(dim() == other.dim());

   // If the result has not been allocated then do so.  Otherwise size it
   // correctly.
   if (result == 0) {
      result = new Vector(d_dim, d_distributed);
   }
   else {
      result->setSize(d_dim);
   }

   // Do the addition.
   for (int i = 0; i < d_dim; ++i) {
      result->d_vec[i] = d_vec[i] + factor*other.d_vec[i];
   }
}

void
Vector::plusAx(
   double factor,
   const Vector& other,
   Vector& result) const
{
   CAROM_ASSERT(result.distributed() == distributed());
   CAROM_ASSERT(distributed() == other.distributed());
   CAROM_ASSERT(dim() == other.dim());

   // Size result correctly.
   result.setSize(d_dim);

   // Do the addition.
   for (int i = 0; i < d_dim; ++i) {
      result.d_vec[i] = d_vec[i] + factor*other.d_vec[i];
   }
}

void
Vector::plusEqAx(
   double factor,
   const Vector& other)
{
   CAROM_ASSERT(distributed() == other.distributed());
   CAROM_ASSERT(dim() == other.dim());

   // Do the addition.
   for (int i = 0; i < d_dim; ++i) {
      d_vec[i] += factor*other.d_vec[i];
   }
}

void
Vector::minus(
   const Vector& other,
   Vector*& result) const
{
   CAROM_ASSERT(result == 0 || result->distributed() == distributed());
   CAROM_ASSERT(distributed() == other.distributed());
   CAROM_ASSERT(dim() == other.dim());

   // If the result has not been allocated then do so.  Otherwise size it
   // correctly.
   if (result == 0) {
      result = new Vector(d_dim, d_distributed);
   }
   else {
      result->setSize(d_dim);
   }

   // Do the subtraction.
   for (int i = 0; i < d_dim; ++i) {
      result->d_vec[i] = d_vec[i] - other.d_vec[i];
   }
}

void
Vector::minus(
   const Vector& other,
   Vector& result) const
{
   CAROM_ASSERT(result.distributed() == distributed());
   CAROM_ASSERT(distributed() == other.distributed());
   CAROM_ASSERT(dim() == other.dim());

   // Size result correctly.
   result.setSize(d_dim);

   // Do the subtraction.
   for (int i = 0; i < d_dim; ++i) {
      result.d_vec[i] = d_vec[i] - other.d_vec[i];
   }
}

void
Vector::mult(
   double factor,
   Vector*& result) const
{
   CAROM_ASSERT(result == 0 || result->distributed() == distributed());

   // If the result has not been allocated then do so.  Otherwise size it
   // correctly.
   if (result == 0) {
      result = new Vector(d_dim, d_distributed);
   }
   else {
      result->setSize(d_dim);
   }

   // Do the multiplication.
   for (int i = 0; i < d_dim; ++i) {
      result->d_vec[i] = factor*d_vec[i];
   }
}

void
Vector::mult(
   double factor,
   Vector& result) const
{
   CAROM_ASSERT(result.distributed() == distributed());

   // Size result correctly.
   result.setSize(d_dim);

   // Do the multiplication.
   for (int i = 0; i < d_dim; ++i) {
      result.d_vec[i] = factor*d_vec[i];
   }
}

void
Vector::write(const std::string& base_file_name)
{
   CAROM_ASSERT(!base_file_name.empty());

   int mpi_init;
   MPI_Initialized(&mpi_init);
   int rank;
   if (mpi_init) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   }
   else {
      rank = 0;
   }

   char tmp[100];
   sprintf(tmp, ".%06d", rank);
   std::string full_file_name = base_file_name + tmp;
   HDFDatabase database;
   database.create(full_file_name);

   sprintf(tmp, "distributed");
   database.putInteger(tmp, d_distributed);
   sprintf(tmp, "dim");
   database.putInteger(tmp, d_dim);
   sprintf(tmp, "data");
   database.putDoubleArray(tmp, d_vec, d_dim);
   database.close();
}

void
Vector::print(const char * prefix)
{
   int my_rank;
   const bool success = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   CAROM_ASSERT(success);

   std::string filename_str = prefix + std::to_string(my_rank);
   const char * filename = filename_str.c_str();
   FILE * pFile = fopen(filename,"w");
   for (int k = 0; k < d_dim; ++k) {
     fprintf(pFile, " %25.20e\n", d_vec[k]);
   }
   fclose(pFile);
}

void
Vector::read(const std::string& base_file_name)
{
   CAROM_ASSERT(!base_file_name.empty());

   int mpi_init;
   MPI_Initialized(&mpi_init);
   int rank;
   if (mpi_init) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   }
   else {
      rank = 0;
   }

   char tmp[100];
   sprintf(tmp, ".%06d", rank);
   std::string full_file_name = base_file_name + tmp;
   HDFDatabase database;
   database.open(full_file_name);

   sprintf(tmp, "distributed");
   int distributed;
   database.getInteger(tmp, distributed);
   d_distributed = bool(distributed);
   int dim;
   sprintf(tmp, "dim");
   database.getInteger(tmp, dim);
   setSize(dim);
   sprintf(tmp, "data");
   database.getDoubleArray(tmp, d_vec, d_alloc_size);
   d_owns_data = true;
   if (mpi_init) {
      MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
   }
   else {
      d_num_procs = 1;
   }
   database.close();
}

}
