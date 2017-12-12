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

// Description: A simple, parallel Matrix class with the utility needed to
//              support the basis generation methods of this library.  A
//              distributed Matrix has its rows distributed across processors.

#include "Matrix.h"

#include "mpi.h"
#include <string.h>

extern "C" {
// LU decomposition of a general matrix.
void
dgetrf_(int*, int*, double*, int*, int*, int*);

// Generate inverse of a matrix given its LU decomposition.
void
dgetri_(int*, double*, int*, int*, double*, int*, int*);
}

namespace CAROM {

Matrix::Matrix(
   int num_rows,
   int num_cols,
   bool distributed) :
   d_mat(0),
   d_alloc_size(0),
   d_distributed(distributed)
{
   CAROM_ASSERT(num_rows > 0);
   CAROM_ASSERT(num_cols > 0);
   setSize(num_rows, num_cols);
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init) {
      MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
   }
   else {
      d_num_procs = 1;
   }
}

Matrix::Matrix(
   const double* mat,
   int num_rows,
   int num_cols,
   bool distributed) :
   d_mat(0),
   d_alloc_size(0),
   d_distributed(distributed)
{
   CAROM_ASSERT(mat != 0);
   CAROM_ASSERT(num_rows > 0);
   CAROM_ASSERT(num_cols > 0);
   setSize(num_rows, num_cols);
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init) {
      MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
   }
   else {
      d_num_procs = 1;
   }
   memcpy(d_mat, mat, d_alloc_size*sizeof(double));
}

Matrix::Matrix(
   const Matrix& other) :
   d_mat(0),
   d_alloc_size(0),
   d_distributed(other.d_distributed)
{
   setSize(other.d_num_rows, other.d_num_cols);
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init) {
      MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
   }
   else {
      d_num_procs = 1;
   }
   memcpy(d_mat, other.d_mat, d_alloc_size*sizeof(double));
}

Matrix::~Matrix()
{
   if (d_mat) {
      delete [] d_mat;
   }
}

Matrix&
Matrix::operator = (
   const Matrix& rhs)
{
   d_distributed = rhs.d_distributed;
   d_num_procs = rhs.d_num_procs;
   setSize(rhs.d_num_rows, rhs.d_num_cols);
   memcpy(d_mat, rhs.d_mat, d_num_rows*d_num_cols*sizeof(double));
   return *this;
}

void
Matrix::mult(
   const Matrix& other,
   Matrix*& result) const
{
   CAROM_ASSERT(result == 0 || result->distributed() == distributed());
   CAROM_ASSERT(!other.distributed());
   CAROM_ASSERT(numColumns() == other.numRows());

   // If the result has not been allocated then do so.  Otherwise size it
   // correctly.
   if (result == 0) {
      result = new Matrix(d_num_rows, other.d_num_cols, d_distributed);
   }
   else {
      result->setSize(d_num_rows, other.d_num_cols);
   }
   for (int this_row = 0; this_row < d_num_rows; ++this_row) {
      for (int other_col = 0; other_col < other.d_num_cols; ++other_col) {
         double result_val = 0.0;
         for (int entry = 0; entry < d_num_cols; ++entry) {
            result_val += item(this_row, entry)*other.item(entry, other_col);
         }
         result->item(this_row, other_col) = result_val;
      }
   }
}

void
Matrix::mult(
   const Vector& other,
   Vector*& result) const
{
   CAROM_ASSERT(result == 0 || result->distributed() == distributed());
   CAROM_ASSERT(!other.distributed());
   CAROM_ASSERT(numColumns() == other.dim());

   // If the result has not been allocated then do so.  Otherwise size it
   // correctly.
   if (result == 0) {
      result = new Vector(d_num_rows, d_distributed);
   }
   else {
      result->setSize(d_num_rows);
   }
   for (int this_row = 0; this_row < d_num_rows; ++this_row) {
      double result_val = 0.0;
      for (int entry = 0; entry < d_num_cols; ++entry) {
         result_val += item(this_row, entry)*other.item(entry);
      }
      result->item(this_row) = result_val;
   }
}

void
Matrix::multPlus(
   Vector& a,
   const Vector& b,
   double c) const
{
   CAROM_ASSERT(a.distributed() == distributed());
   CAROM_ASSERT(!b.distributed());
   CAROM_ASSERT(numColumns() == b.dim());
   CAROM_ASSERT(numRows() == a.dim());

   for (int this_row = 0; this_row < d_num_rows; ++this_row) {
      double tmp = 0.0;
      for (int this_col = 0; this_col < d_num_cols; ++this_col) {
         tmp += item(this_row, this_col)*b.item(this_col);
      }
      a.item(this_row) += tmp*c;
   }
}

void
Matrix::transposeMult(
   const Matrix& other,
   Matrix*& result) const
{
   CAROM_ASSERT(result == 0 || !result->distributed());
   CAROM_ASSERT(distributed() == other.distributed());
   CAROM_ASSERT(numRows() == other.numRows());

   // If the result has not been allocated then do so.  Otherwise size it
   // correctly.
   if (result == 0) {
      result = new Matrix(d_num_cols, other.d_num_cols, false);
   }
   else {
      result->setSize(d_num_cols, other.d_num_cols);
   }
   for (int this_col = 0; this_col < d_num_cols; ++this_col) {
      for (int other_col = 0; other_col < other.d_num_cols; ++other_col) {
         double result_val = 0.0;
         for (int entry = 0; entry < d_num_rows; ++entry) {
            result_val += item(entry, this_col)*other.item(entry, other_col);
         }
         result->item(this_col, other_col) = result_val;
      }
   }
   if (d_distributed && d_num_procs > 1) {
      int new_mat_size = d_num_cols*other.d_num_cols;
      MPI_Allreduce(MPI_IN_PLACE,
         &result->item(0, 0),
         new_mat_size,
         MPI_DOUBLE,
         MPI_SUM,
         MPI_COMM_WORLD);
   }
}

void
Matrix::transposeMult(
   const Vector& other,
   Vector*& result) const
{
   CAROM_ASSERT(result == 0 || !result->distributed());
   CAROM_ASSERT(distributed() == other.distributed());
   CAROM_ASSERT(numRows() == other.dim());

   // If the result has not been allocated then do so.  Otherwise size it
   // correctly.
   if (result == 0) {
      result = new Vector(d_num_cols, false);
   }
   else {
      result->setSize(d_num_cols);
   }
   for (int this_col = 0; this_col < d_num_cols; ++this_col) {
      double result_val = 0.0;
      for (int entry = 0; entry < d_num_rows; ++entry) {
         result_val += item(entry, this_col)*other.item(entry);
      }
      result->item(this_col) = result_val;
   }
   if (d_distributed && d_num_procs > 1) {
      MPI_Allreduce(MPI_IN_PLACE,
         &result->item(0),
         d_num_cols,
         MPI_DOUBLE,
         MPI_SUM,
         MPI_COMM_WORLD);
   }
}

void
Matrix::inverse(
   Matrix*& result) const
{
   CAROM_ASSERT(result == 0 ||
                (!result->distributed() &&
                 result->numRows() == numRows() &&
                 result->numColumns() == numColumns()));
   CAROM_ASSERT(!distributed());
   CAROM_ASSERT(numRows() == numColumns());

   // If the result has not been allocated then do so.  Otherwise size it
   // correctly.
   if (result == 0) {
      result = new Matrix(d_num_rows, d_num_cols, false);
   }
   else {
      result->setSize(d_num_rows, d_num_cols);
   }

   // Call lapack routines to do the inversion.
   // Set up some stuff the lapack routines need.
   int info;
   int mtx_size = d_num_rows;
   int lwork = mtx_size*mtx_size;
   int* ipiv = new int [mtx_size];
   double* work = new double [lwork];
   // To use lapack we need a column major representation of this which is
   // essentially the transform of this.  Use result for this representation.
   for (int row = 0; row < mtx_size; ++row) {
      for (int col = 0; col < mtx_size; ++col) {
         result->item(col, row) = item(row, col);
      }
   }
   // Now call lapack to do the inversion.
   dgetrf_(&mtx_size, &mtx_size, result->d_mat, &mtx_size, ipiv, &info);
   dgetri_(&mtx_size, result->d_mat, &mtx_size, ipiv, work, &lwork, &info);
   // Result now has the inverse in a column major representation.  Put it
   // into row major order.
   for (int row = 0; row < mtx_size; ++row) {
      for (int col = row+1; col < mtx_size; ++col) {
         double tmp = result->item(row, col);
         result->item(row, col) = result->item(col, row);
         result->item(col, row) = tmp;
      }
   }
}

void
Matrix::inverse()
{
   CAROM_ASSERT(!distributed());
   CAROM_ASSERT(numRows() == numColumns());

   // Call lapack routines to do the inversion.
   // Set up some stuff the lapack routines need.
   int info;
   int mtx_size = d_num_rows;
   int lwork = mtx_size*mtx_size;
   int* ipiv = new int [mtx_size];
   double* work = new double [lwork];
   // To use lapack we need a column major representation of this which is
   // essentially the transform of this.
   for (int row = 0; row < mtx_size; ++row) {
      for (int col = row+1; col < mtx_size; ++col) {
         double tmp = item(row, col);
         item(row, col) = item(col, row);
         item(col, row) = tmp;
      }
   }
   // Now call lapack to do the inversion.
   dgetrf_(&mtx_size, &mtx_size, d_mat, &mtx_size, ipiv, &info);
   dgetri_(&mtx_size, d_mat, &mtx_size, ipiv, work, &lwork, &info);
   // This now has its inverse in a column major representation.  Put it into
   // row major representation.
   for (int row = 0; row < mtx_size; ++row) {
      for (int col = row+1; col < mtx_size; ++col) {
         double tmp = item(row, col);
         item(row, col) = item(col, row);
         item(col, row) = tmp;
      }
   }
}

}
