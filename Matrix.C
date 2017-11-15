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

namespace CAROM {

Matrix::Matrix(
   int num_rows,
   int num_cols,
   bool distributed) :
   d_num_rows(num_rows),
   d_num_cols(num_cols),
   d_distributed(distributed)
{
   CAROM_ASSERT(num_rows > 0);
   CAROM_ASSERT(num_cols > 0);
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init) {
      MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
   }
   else {
      d_num_procs = 1;
   }
   d_mat = new double [num_cols*num_rows];
}

Matrix::Matrix(
   const double* mat,
   int num_rows,
   int num_cols,
   bool distributed) :
   d_num_rows(num_rows),
   d_num_cols(num_cols),
   d_distributed(distributed)
{
   CAROM_ASSERT(mat != 0);
   CAROM_ASSERT(num_rows > 0);
   CAROM_ASSERT(num_cols > 0);
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init) {
      MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
   }
   else {
      d_num_procs = 1;
   }
   int num_entries = d_num_cols*d_num_rows;
   d_mat = new double [num_entries];
   memcpy(d_mat, mat, num_entries*sizeof(double));
}

Matrix::Matrix(
   const Matrix& other) :
   d_num_rows(other.d_num_rows),
   d_num_cols(other.d_num_cols),
   d_distributed(other.d_distributed)
{
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init) {
      MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
   }
   else {
      d_num_procs = 1;
   }
   int num_entries = d_num_cols*d_num_rows;
   d_mat = new double [num_entries];
   memcpy(d_mat, other.d_mat, num_entries*sizeof(double));
}

Matrix::~Matrix()
{
   delete [] d_mat;
}

Matrix&
Matrix::operator = (
   const Matrix& rhs)
{
   bool realloc = (d_num_cols*d_num_rows) != (rhs.d_num_cols*rhs.d_num_rows);
   d_num_rows = rhs.d_num_rows;
   d_num_cols = rhs.d_num_cols;
   d_distributed = rhs.d_distributed;
   if (realloc) {
      delete [] d_mat;
      d_mat = new double[d_num_cols*d_num_rows];
   }
   memcpy(d_mat, rhs.d_mat, d_num_cols*d_num_rows*sizeof(double));
   return *this;
}

Matrix*
Matrix::mult(
   const Matrix& other) const
{
   CAROM_ASSERT(!other.distributed());
   CAROM_ASSERT(numColumns() == other.numRows());
   Matrix* result = new Matrix(d_num_rows,
                               other.d_num_cols,
                               d_distributed);
   for (int this_row = 0; this_row < d_num_rows; ++this_row) {
      for (int other_col = 0; other_col < other.d_num_cols; ++other_col) {
         double result_val = 0.0;
         for (int entry = 0; entry < d_num_cols; ++entry) {
            result_val += item(this_row, entry)*other.item(entry, other_col);
         }
         result->item(this_row, other_col) = result_val;
      }
   }
   return result;
}

Vector*
Matrix::mult(
   const Vector& other) const
{
   CAROM_ASSERT(distributed() && !other.distributed());
   CAROM_ASSERT(numColumns() == other.dim());
   Vector* result = new Vector(d_num_rows, true);
   for (int this_row = 0; this_row < d_num_rows; ++this_row) {
      double result_val = 0.0;
      for (int entry = 0; entry < d_num_cols; ++entry) {
         result_val += item(this_row, entry)*other.item(entry);
      }
      result->item(this_row) = result_val;
   }
   return result;
}

Matrix*
Matrix::transposeMult(
   const Matrix& other) const
{
   CAROM_ASSERT(distributed() == other.distributed());
   CAROM_ASSERT(numRows() == other.numRows());
   Matrix* result = new Matrix(d_num_cols, other.d_num_cols, false);
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
   return result;
}

Vector*
Matrix::transposeMult(
   const Vector& other) const
{
   CAROM_ASSERT(distributed() && other.distributed());
   CAROM_ASSERT(numRows() == other.dim());
   Vector* result = new Vector(d_num_cols, false);
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
   return result;
}

}
