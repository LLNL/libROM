/******************************************************************************
 *
 * Copyright (c) 2013-2016, Lawrence Livermore National Security, LLC.
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
