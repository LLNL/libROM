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

// BLAS-3 version of QR decomposition with column pivoting
void
dgeqp3_(int*, int*, double*, int*, int*, double*, double*, int*, int*);
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

bool
Matrix::balanced() const
{
  // TODO(oxberry1@llnl.gov): Relax the assumption in libROM that
  // objects use MPI_COMM_WORLD, and that rank 0 is the master rank of
  // an object.
  const int master_rank = 0;
  const MPI_Comm comm = MPI_COMM_WORLD;

  // A Matrix is "balanced" (load-balanced for distributed dense matrix
  // computations) if:
  //
  // (1) the number of rows owned by each process in any pair of
  // processes on the communicator differs by at most one
  //
  // (2) process j has no fewer rows than k if j is less than k (j and
  // k are both integers corresponding to process ranks)

  // Serial matrices are balanced by definition; one process owns all
  // rows
  if (!distributed()) return true;

  // Otherwise, get the total number of rows of the matrix.
  int num_total_rows = d_num_rows;
  const int reduce_count = 1;
  CAROM_ASSERT(MPI_Allreduce(MPI_IN_PLACE,
			     &num_total_rows,
			     reduce_count,
			     MPI_INT,
			     MPI_SUM,
			     comm) == MPI_SUCCESS);

  const int first_rank_with_fewer = num_total_rows % d_num_procs;
  int my_rank;
  CAROM_ASSERT(MPI_Comm_rank(comm, &my_rank) == MPI_SUCCESS);

  const int  min_rows_per_rank = num_total_rows / d_num_procs;
  const bool has_extra_row     = my_rank < first_rank_with_fewer;
  const int  max_rows_on_rank  = min_rows_per_rank + has_extra_row;
  const bool has_enough_rows   = (d_num_rows >= min_rows_per_rank);
  const bool has_too_many_rows = (d_num_rows > max_rows_on_rank);

  int result = (has_enough_rows && !has_too_many_rows);
  CAROM_ASSERT(MPI_Allreduce(MPI_IN_PLACE,
			     &result,
			     reduce_count,
			     MPI_INT,
			     MPI_LAND,
			     comm) == MPI_SUCCESS);

  return result;
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

void
Matrix::qrcp_pivots_transpose(int* row_pivot,
			      int* row_pivot_owner,
			      int  pivots_requested) const
{
  // For now, implement serial version
  CAROM_ASSERT(!distributed());
  qrcp_pivots_transpose_serial(leading_pivots);
}

void
Matrix::qrcp_pivots_transpose_serial(std::vector<int>& leading_pivots) const
{
  // This method assumes this matrix is serial
  CAROM_ASSERT(!distributed());

  // Number of pivots requested can't exceed the number of rows of the
  // matrix
  CAROM_ASSERT(pivots_requested <= numRows());
  CAROM_ASSERT(pivots_requested > 0);

  // Make sure arrays are allocated before entry; this method does not
  // own the input pointers
  CAROM_ASSERT(row_pivot != NULL);
  CAROM_ASSERT(row_pivot_owner != NULL);

  // Get dimensions of transpose of matrix
  int num_rows_of_transpose = numColumns();
  int num_cols_of_transpose = numRows();

  // LAPACK routines tend to overwrite their inputs, but we'd like to
  // keep the basis matrix and use it in later computations, so copy
  // the basis matrix here.
  Matrix scratch(*this);

  // Allocate work arrays; work array for QR must be at least 1 plus 3
  // times the number of columns of its matrix input. This algorithm
  // applies QR to transposed basis matrix, so the applicable
  // dimension is the number of rows of the basis matrix. It's
  // possible to get better performance by computing the optimal block
  // size and then using that value to size the work array; see the
  // LAPACK source code and documentation for details.
  int lwork = 20 * num_cols_of_transpose + 1;
  double* work = new double[lwork];
  double* tau  = new double[std::min(num_rows_of_transpose,
				     num_cols_of_transpose)];
  int* pivot = new int[num_cols_of_transpose];
  int info;

   // Compute the QR decomposition with column pivots of the transpose
   // of this matrix by abusing the fact that the C++ memory model is
   // row-major format, which is the transpose of the Fortran memory
   // model (which is column-major). Passing the row-major data
   // looks like an in-place transposition to Fortran.
   dgeqp3_(&num_rows_of_transpose,
	   &num_cols_of_transpose,
	   scratch.d_mat,
	   &num_rows_of_transpose,
	   pivot,
	   tau,
	   work,
	   &lwork,
	   &info);

   // Fail if error in LAPACK routine.
   CAROM_ASSERT(info == 0);

   // Assume communicator is MPI_COMM_WORLD and get rank of this
   // process
   int is_mpi_initialized, is_mpi_finalized;
   CAROM_ASSERT(MPI_Initialized(&is_mpi_initialized) == MPI_SUCCESS);
   CAROM_ASSERT(MPI_Finalized(&is_mpi_finalized) == MPI_SUCCESS);
   int my_rank = 0;
   if(is_mpi_initialized && !is_mpi_finalized) {
     const MPI_Comm my_comm = MPI_COMM_WORLD;
     CAROM_ASSERT(MPI_Comm_rank(my_comm, &my_rank) == MPI_SUCCESS);
   }

   // Copy over pivots and subtract one to convert them from a
   // Fortran-based indexing convention (first element of 1-D array by
   // default corresponds to index of 1, though this convention can be
   // overridden) to a C-based indexing convention (first element of
   // 1-D array corresponds to index of 0).
   for (int i = 0; i < pivots_requested; i++) {
     row_pivot[i]       = pivot[i] - 1;
     row_pivot_owner[i] = my_rank;
   }

   // Free arrays
   delete [] work;
   delete [] tau;
   delete [] pivot;
}

}
