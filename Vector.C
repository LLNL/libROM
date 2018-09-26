/******************************************************************************
 *
 * Copyright (c) 2013-2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 * Written by William Arrighi wjarrighi@llnl.gov
 * CODE-686965
 * All rights reserved.
 *
 * This file is part of libROM.
 * For details, see https://computation.llnl.gov/librom
 * Please also read README_BSD_NOTICE.
 *
 * Redistribution and use in source and binary forms, with or without
 * modifications, are permitted provided that the following conditions are met:
 *
 *    o Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the disclaimer below.
 *    o Redistribution in binary form must reproduce the above copyright
 *      notice, this list of conditions and the disclaimer (as noted below) in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *    o Neither the name of the LLNS/LLNL nor the names of its contributors may
 *      be used to endorse or promote products derived from this software
 *      without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
 * LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OR SUCH DAMAGE.
 *
 *****************************************************************************/

// Description: A simple, parallel Vector class with the utility needed to
//              support the basis generation methods of this library.  A
//              distributed Vector has its rows distributed across processors.

#include "Vector.h"

#include "mpi.h"

#include <cmath>
#include <string.h>

namespace CAROM {

Vector::Vector(
   int dim,
   bool distributed) :
   d_vec(0),
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
   d_vec(0),
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
   d_vec(0),
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

void
Vector::zero()
{
   for(int i=0; i<d_dim; ++i) d_vec[i] = 0.0;
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

}
