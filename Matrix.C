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
#include "HDFDatabase.h"

#include "mpi.h"
#include <string.h>

#ifdef CAROM_HAS_ELEMENTAL
#include <El.hpp>
#endif

/* Use Autotools-detected Fortran name-mangling scheme */
#define dgetrf FC_FUNC(dgetrf, DGETRF)
#define dgetri FC_FUNC(dgetri, DGETRI)
#define dgeqp3 FC_FUNC(dgeqp3, DGEQP3)

extern "C" {
// LU decomposition of a general matrix.
void
dgetrf(int*, int*, double*, int*, int*, int*);

// Generate inverse of a matrix given its LU decomposition.
void
dgetri(int*, double*, int*, int*, double*, int*, int*);

// BLAS-3 version of QR decomposition with column pivoting
void
dgeqp3(int*, int*, double*, int*, int*, double*, double*, int*, int*);
}

namespace CAROM {

Matrix::Matrix() :
   d_mat(NULL),
   d_alloc_size(0),
   d_distributed(false),
   d_owns_data(true)
{}

Matrix::Matrix(
   int num_rows,
   int num_cols,
   bool distributed) :
   d_mat(0),
   d_alloc_size(0),
   d_distributed(distributed),
   d_owns_data(true)
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
   double* mat,
   int num_rows,
   int num_cols,
   bool distributed,
   bool copy_data) :
   d_mat(0),
   d_alloc_size(0),
   d_distributed(distributed),
   d_owns_data(copy_data)
{
   CAROM_ASSERT(mat != 0);
   CAROM_ASSERT(num_rows > 0);
   CAROM_ASSERT(num_cols > 0);
   if (copy_data) {
      setSize(num_rows, num_cols);
      memcpy(d_mat, mat, d_alloc_size*sizeof(double));
   }
   else {
      d_mat = mat;
      d_alloc_size = num_rows*num_cols;
      d_num_cols = num_cols;
      d_num_rows = num_rows;
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

Matrix::Matrix(
   const Matrix& other) :
   d_mat(0),
   d_alloc_size(0),
   d_distributed(other.d_distributed),
   d_owns_data(true)
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
   if (d_owns_data && d_mat) {
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

Matrix&
Matrix::operator += (
   const Matrix& rhs)
{
   CAROM_ASSERT(rhs.d_num_rows == d_num_rows);
   CAROM_ASSERT(rhs.d_num_cols == d_num_cols);
   for(int i=0; i<d_num_rows*d_num_cols; ++i) d_mat[i] += rhs.d_mat[i];
   return *this;
}

Matrix&
Matrix::operator -= (
   const Matrix& rhs)
{
   CAROM_ASSERT(rhs.d_num_rows == d_num_rows);
   CAROM_ASSERT(rhs.d_num_cols == d_num_cols);
   for(int i=0; i<d_num_rows*d_num_cols; ++i) d_mat[i] -= rhs.d_mat[i];
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

Matrix&
Matrix::operator = (
   const double a)
{
   for(int i=0; i<d_num_rows*d_num_cols; ++i) {
     d_mat[i] = a;
   }
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

   // Do the multiplication.
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
   const Matrix& other,
   Matrix& result) const
{
   CAROM_ASSERT(result.distributed() == distributed());
   CAROM_ASSERT(!other.distributed());
   CAROM_ASSERT(numColumns() == other.numRows());

   // Size result correctly.
   result.setSize(d_num_rows, other.d_num_cols);

   // Do the multiplication.
   for (int this_row = 0; this_row < d_num_rows; ++this_row) {
      for (int other_col = 0; other_col < other.d_num_cols; ++other_col) {
         double result_val = 0.0;
         for (int entry = 0; entry < d_num_cols; ++entry) {
            result_val += item(this_row, entry)*other.item(entry, other_col);
         }
         result.item(this_row, other_col) = result_val;
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

   // Do the multiplication.
   for (int this_row = 0; this_row < d_num_rows; ++this_row) {
      double result_val = 0.0;
      for (int entry = 0; entry < d_num_cols; ++entry) {
         result_val += item(this_row, entry)*other.item(entry);
      }
      result->item(this_row) = result_val;
   }
}

void
Matrix::mult(
   const Vector& other,
   Vector& result) const
{
   CAROM_ASSERT(result.distributed() == distributed());
   CAROM_ASSERT(!other.distributed());
   CAROM_ASSERT(numColumns() == other.dim());

   // Size result correctly.
   result.setSize(d_num_rows);

   // Do the multiplication.
   for (int this_row = 0; this_row < d_num_rows; ++this_row) {
      double result_val = 0.0;
      for (int entry = 0; entry < d_num_cols; ++entry) {
         result_val += item(this_row, entry)*other.item(entry);
      }
      result.item(this_row) = result_val;
   }
}

void
Matrix::pointwise_mult(
   int this_row,
   const Vector& other,
   Vector& result) const
{
   CAROM_ASSERT(!result.distributed());
   CAROM_ASSERT(!distributed());
   CAROM_ASSERT(!other.distributed());
   CAROM_ASSERT(numColumns() == other.dim());
 
   // Do the multiplication.
   for (int entry = 0; entry < d_num_cols; ++entry) {
     result.item(entry) = item(this_row, entry)*other.item(entry);
   }
}

void
Matrix::pointwise_mult(
   int this_row,
   Vector& other) const
{
   CAROM_ASSERT(!distributed());
   CAROM_ASSERT(!other.distributed());
   CAROM_ASSERT(numColumns() == other.dim());
 
   // Do the multiplication.
   for (int entry = 0; entry < d_num_cols; ++entry) {
     other.item(entry) *= item(this_row, entry);
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

   // Do the multiplication.
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
   const Matrix& other,
   Matrix& result) const
{
   CAROM_ASSERT(!result.distributed());
   CAROM_ASSERT(distributed() == other.distributed());
   CAROM_ASSERT(numRows() == other.numRows());

   // Size result correctly.
   result.setSize(d_num_cols, other.d_num_cols);

   // Do the multiplication.
   for (int this_col = 0; this_col < d_num_cols; ++this_col) {
      for (int other_col = 0; other_col < other.d_num_cols; ++other_col) {
         double result_val = 0.0;
         for (int entry = 0; entry < d_num_rows; ++entry) {
            result_val += item(entry, this_col)*other.item(entry, other_col);
         }
         result.item(this_col, other_col) = result_val;
      }
   }
   if (d_distributed && d_num_procs > 1) {
      int new_mat_size = d_num_cols*other.d_num_cols;
      MPI_Allreduce(MPI_IN_PLACE,
         &result.item(0, 0),
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

   // Do the multiplication.
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
Matrix::transposeMult(
   const Vector& other,
   Vector& result) const
{
   CAROM_ASSERT(!result.distributed());
   CAROM_ASSERT(distributed() == other.distributed());
   CAROM_ASSERT(numRows() == other.dim());

   // If the result has not been allocated then do so.  Otherwise size it
   // correctly.
   result.setSize(d_num_cols);

   // Do the multiplication.
   for (int this_col = 0; this_col < d_num_cols; ++this_col) {
      double result_val = 0.0;
      for (int entry = 0; entry < d_num_rows; ++entry) {
         result_val += item(entry, this_col)*other.item(entry);
      }
      result.item(this_col) = result_val;
   }
   if (d_distributed && d_num_procs > 1) {
      MPI_Allreduce(MPI_IN_PLACE,
         &result.item(0),
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
   dgetrf(&mtx_size, &mtx_size, result->d_mat, &mtx_size, ipiv, &info);
   dgetri(&mtx_size, result->d_mat, &mtx_size, ipiv, work, &lwork, &info);
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
Matrix::inverse(
   Matrix& result) const
{
   CAROM_ASSERT(!result.distributed() && result.numRows() == numRows() &&
                result.numColumns() == numColumns());
   CAROM_ASSERT(!distributed());
   CAROM_ASSERT(numRows() == numColumns());

   // Size result correctly.
   result.setSize(d_num_rows, d_num_cols);

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
         result.item(col, row) = item(row, col);
      }
   }
   // Now call lapack to do the inversion.
   dgetrf(&mtx_size, &mtx_size, result.d_mat, &mtx_size, ipiv, &info);
   dgetri(&mtx_size, result.d_mat, &mtx_size, ipiv, work, &lwork, &info);
   // Result now has the inverse in a column major representation.  Put it
   // into row major order.
   for (int row = 0; row < mtx_size; ++row) {
      for (int col = row+1; col < mtx_size; ++col) {
         double tmp = result.item(row, col);
         result.item(row, col) = result.item(col, row);
         result.item(col, row) = tmp;
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
   dgetrf(&mtx_size, &mtx_size, d_mat, &mtx_size, ipiv, &info);
   dgetri(&mtx_size, d_mat, &mtx_size, ipiv, work, &lwork, &info);
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

void Matrix::transposePseudoinverse()
{
   CAROM_ASSERT(!distributed());
   CAROM_ASSERT(numRows() >= numColumns());

   if (numRows() == numColumns())
     {
       inverse();
     }
   else
     {
       Matrix *AtA = this->transposeMult(this);

       // Directly invert AtA, which is a bad idea if AtA is not small. 
       AtA->inverse();

       // Pseudoinverse is (AtA)^{-1}*this^T, but we store the transpose of the result in this, namely this*(AtA)^{-T}.
       Vector row(numColumns(), false);
       Vector res(numColumns(), false);
       for (int i=0; i<numRows(); ++i)
	 { // Compute i-th row of this multiplied by (AtA)^{-T}, whose transpose is (AtA)^{-1} times i-th row transposed.
	   for (int j=0; j<numColumns(); ++j)
	     row.item(j) = this->item(i,j);
	   
	   AtA->mult(row, res);

	   // Overwrite i-th row with transpose of result.
	   for (int j=0; j<numColumns(); ++j)
	     this->item(i,j) = res.item(j);
	 }
   
       delete AtA;
     }
}

void
Matrix::print(const char * prefix)
{
   int my_rank;
   const bool success = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   CAROM_ASSERT(success);

   std::string filename_str = prefix + std::to_string(my_rank); 
   const char * filename = filename_str.c_str();
   FILE * pFile = fopen(filename,"w");
   for (int row = 0; row < d_num_rows; ++row) { 
     for (int col = 0; col < d_num_cols; ++col) { 
       fprintf(pFile, " %25.20e\t", item(row,col));
     }
     fprintf(pFile, "\n");
   }
   fclose(pFile);
}

void
Matrix::write(const std::string& base_file_name)
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
   sprintf(tmp, "num_rows");
   database.putInteger(tmp, d_num_rows);
   sprintf(tmp, "num_cols");
   database.putInteger(tmp, d_num_cols);
   sprintf(tmp, "data");
   database.putDoubleArray(tmp, d_mat, d_num_rows*d_num_cols);
   database.close();
}

void
Matrix::read(const std::string& base_file_name)
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
   int num_rows;
   sprintf(tmp, "num_rows");
   database.getInteger(tmp, num_rows);
   int num_cols;
   sprintf(tmp, "num_cols");
   database.getInteger(tmp, num_cols);
   setSize(num_rows,num_cols);
   sprintf(tmp, "data");
   database.getDoubleArray(tmp, d_mat, d_alloc_size);
   d_owns_data = true;
   if (mpi_init) {
      MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
   }
   else {
      d_num_procs = 1;
   }
   database.close();
}

void
Matrix::qrcp_pivots_transpose(int* row_pivot,
			      int* row_pivot_owner,
			      int  pivots_requested) const
{
  if(!distributed()) {
    return qrcp_pivots_transpose_serial(row_pivot,
					row_pivot_owner,
					pivots_requested);
  }
  else{
#ifdef CAROM_HAS_ELEMENTAL
    return qrcp_pivots_transpose_distributed(row_pivot,
					     row_pivot_owner,
					     pivots_requested);
#else
    CAROM_ASSERT(false);
#endif
  }
}

void
Matrix::qrcp_pivots_transpose_serial(int* row_pivot,
				     int* row_pivot_owner,
				     int  pivots_requested) const
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
   dgeqp3(&num_rows_of_transpose,
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

void
Matrix::qrcp_pivots_transpose_distributed(int* row_pivot,
					  int* row_pivot_owner,
					  int  pivots_requested)
const
{
#ifdef CAROM_HAS_ELEMENTAL
  // Shim to design interface; not implemented yet

  // Check if distributed; otherwise, use serial implementation
  CAROM_ASSERT(distributed());

  // Elemental implementation
  return qrcp_pivots_transpose_distributed_elemental
    (row_pivot, row_pivot_owner, pivots_requested);

  // TODO(oxberry1): ScaLAPACK implementation?
#else
  CAROM_ASSERT(false);
#endif
}

void
Matrix::qrcp_pivots_transpose_distributed_elemental
(int* row_pivot, int* row_pivot_owner, int pivots_requested)
const
{
#ifdef CAROM_HAS_ELEMENTAL
  // Check if distributed; otherwise, use serial implementation
  CAROM_ASSERT(distributed());

  // Check if balanced
  if (balanced()) {
    qrcp_pivots_transpose_distributed_elemental_balanced
      (row_pivot, row_pivot_owner, pivots_requested);
  }
  else {
    qrcp_pivots_transpose_distributed_elemental_unbalanced
      (row_pivot, row_pivot_owner, pivots_requested);
  }
#else
  CAROM_ASSERT(false);
#endif
}

void
Matrix::qrcp_pivots_transpose_distributed_elemental_balanced
(int* row_pivot, int* row_pivot_owner, int pivots_requested)
const
{
#ifdef CAROM_HAS_ELEMENTAL
  // Compute pivots redundantly across all processes using the QRCP
  // from the distributed dense linear algebra library Elemental.

  // The following assumptions are made in this implementation just to
  // get a version up and running:
  //
  // (1) this->balanced() == true; // The matrix is balanced
  //
  // (2) Process 0 is the master rank of the object
  //
  // (3) Matrix rows are distributed block-cyclically over all processes
  //     starting on process zero with row 0, in increasing order of row
  //     index (all elements of a given row are still stored on the same
  //     process)
  //
  // (4) This Matrix is distributed over the MPI_COMM_WORLD
  //     communicator
  //
  // Some of these assumptions can be relaxed if the Matrix object
  // stores more information.

  // Check if distributed and balanced
  CAROM_ASSERT(distributed() && balanced());

  // Make sure arrays are allocated before entry; this method does not
  // own the input pointers
  CAROM_ASSERT(row_pivot != NULL);
  CAROM_ASSERT(row_pivot_owner != NULL);

  // Compute total number of rows to set global sizes of matrix
  const MPI_Comm comm    = MPI_COMM_WORLD;
  const int master_rank  = 0;

  int num_total_rows     = d_num_rows;
  const int reduce_count = 1;
  CAROM_ASSERT(MPI_Allreduce(MPI_IN_PLACE,
			     &num_total_rows,
			     reduce_count,
			     MPI_INT,
			     MPI_SUM,
			     comm) == MPI_SUCCESS);

  // Number of pivots requested can't exceed the total number of rows
  // of the matrix
  CAROM_ASSERT(pivots_requested <= num_total_rows);
  CAROM_ASSERT(pivots_requested > 0);

  // To compute a column-pivoted QRCP using Elemental, we need to
  // first get the data into datatypes Elemental can operate on.

  // Construct a process grid on the communicator that is 1
  // by (# processes), in column-major order; each process in the grid
  // will own its own rows of the matrix
  const El::Grid grid(comm, 1);

  // The row communicator in the grid should have the same number
  // of processes as the comm owned by the Matrix
  CAROM_ASSERT(El::mpi::Size(grid.RowComm()) == d_num_procs);

  // Instantiate the transposed matrix; Elemental calls the number of
  // rows the "height" of the matrix, and the number of columns the
  // "width" of the matrix
  El::Int height = static_cast<El::Int>(numColumns());
  El::Int width  = static_cast<El::Int>(numRows());
  El::Int root   = static_cast<El::Int>(master_rank);
  El::DistMatrix<double> scratch(height, width, grid, root);

  // Columns of the matrix should be distributed round-robin; each
  // element of a column should be on the same process. The
  // distribution should satisfy two invariants:
  //
  // (1) Each "global row" should have the same rank in its column
  // communicator
  //
  // (2) For the scratch matrix global column index j, the process (on
  // the row communicator) that owns j should be process rank (j %
  // d_num_procs).
  //
  CAROM_ASSERT(scratch.RowOwner(0) == scratch.RowOwner(1));
  int my_rank;
  CAROM_ASSERT(MPI_Comm_rank(comm, &my_rank) == MPI_SUCCESS);
  El::Int rank_as_row = static_cast<El::Int>(my_rank);
  CAROM_ASSERT(scratch.ColOwner(rank_as_row) == rank_as_row);

  // Set up work matrices
  El::DistMatrix<double> householder_scalars(grid);
  El::DistMatrix<double> diagonal(grid);
  El::DistPermutation perm(grid);
  El::QRCtrl<double> ctrl;

  // Build the copy of the matrix element by element, mapping
  // local indices to global indices. The Elemental::DistMatrix indices
  // are the transpose of the Matrix indices

  //  Elemental assigns elements of a matrix in a block-cyclic fashion
  // based on the size (specifically, the height and width) of the
  // El::Grid object; see
  // http://libelemental.org/documentation/dev/core/dist_matrix/Element/MC_MR.html
  // for details on the default element distribution, which is what
  // this code uses.
  //
  // If the matrix is balanced, then the matrix rows (i.e., columns of
  // the scratch, transpose matrix) should conform to the row
  // distribution of the Elemental process grid, and additional
  // communication (beyond possibly in local to global index mapping
  // queries) is not needed to assign elements in the proper places
  // because the matrix is already load-balanced. Furthermore, all of
  // the elements of a row of this matrix (column of scratch,
  // transpose matrix) are on a single process. Since the matrix rows
  // (columns of scratch) don't have to be redistributed, we can use
  // Elemental's native (global <--> local) indexing to construct a
  // map between local indices on a given process and global indices
  // of the distributed matrix, and Elemental will also give us the
  // correct process rank owning any given column of the scratch
  // matrix.
  for (int row = 0; row < d_num_rows; row++) {
    El::Int el_loc_col = static_cast<El::Int>(row);
    El::Int el_global_col = scratch.GlobalCol(el_loc_col);
    for (int col = 0; col < d_num_cols; col++) {
      El::Int el_loc_row = static_cast<El::Int>(col);
      El::Int el_global_row = scratch.GlobalRow(el_loc_row);
      scratch.Set(el_global_row, el_global_col, this->item(row, col));
    }
  }

  // After transferring the data over to Elemental's native data
  // types, we compute the QRCP.
  El::QR(scratch, householder_scalars, diagonal, perm); // add ctrl if needed

  // Then, we transfer the pivots into the pivot array
  // stored redundantly on each process.
  for (size_t i = 0; i < pivots_requested; i++) {
    El::Int el_i = static_cast<El::Int>(i);
    El::Int el_perm_i = perm.Image(el_i);

    // The permutation is computed in terms of global indices, so
    // we need to compute the local column pivot of the transpose;
    // this will be a row pivot
    El::Int el_loc_i = scratch.LocalCol(el_perm_i);
    int loc_i = static_cast<int>(el_loc_i);
    row_pivot[i] = loc_i;

    // The global index of the permutation can also be used to figure
    // out which process owns that pivot row, because this process is
    // also the process that owns that global column of the scratch
    // matrix
    El::Int el_owner = scratch.ColOwner(el_perm_i);
    int owner = static_cast<int>(el_owner);
    row_pivot_owner[i] = owner;
  }
#else
  CAROM_ASSERT(false);
#endif
}

void
Matrix::qrcp_pivots_transpose_distributed_elemental_unbalanced
(int* row_pivot, int* row_pivot_owner, int pivots_requested)
const
{
#ifdef CAROM_HAS_ELEMENTAL
  // Compute pivots redundantly across all processes using the QRCP
  // from the distributed dense linear algebra library Elemental.

  // The following assumptions are made in this implementation just to
  // get a version up and running:
  //
  // (1) Process 0 is the master rank of the object
  //
  // (2) This Matrix is distributed over the MPI_COMM_WORLD
  //     communicator
  //
  // Some of these assumptions can be relaxed if the Matrix object
  // stores more information.

  // Check if distributed and unbalanced
  CAROM_ASSERT(distributed() && !balanced());

  // Make sure arrays are allocated before entry; this method does not
  // own the input pointers
  CAROM_ASSERT(row_pivot != NULL);
  CAROM_ASSERT(row_pivot_owner != NULL);

  // Compute total number of rows to set global sizes of matrix
  const MPI_Comm comm    = MPI_COMM_WORLD;
  const int master_rank  = 0;

  int num_total_rows     = d_num_rows;
  const int reduce_count = 1;
  CAROM_ASSERT(MPI_Allreduce(MPI_IN_PLACE,
			     &num_total_rows,
			     reduce_count,
			     MPI_INT,
			     MPI_SUM,
			     comm) == MPI_SUCCESS);

  // Number of pivots requested can't exceed the total number of rows
  // of the matrix
  CAROM_ASSERT(pivots_requested <= num_total_rows);
  CAROM_ASSERT(pivots_requested > 0);

  // To compute a column-pivoted QRCP using Elemental, we need to
  // first get the data into datatypes Elemental can operate on.

  // Construct a process grid on the communicator that is 1
  // by (# processes), in column-major order; each process in the grid
  // will own its own rows of the matrix
  const El::Grid grid(comm, 1);

  // The row communicator in the grid should have the same number
  // of processes as the comm owned by the Matrix
  CAROM_ASSERT(El::mpi::Size(grid.RowComm()) == d_num_procs);

  // Instantiate the transposed matrix; Elemental calls the number of
  // rows the "height" of the matrix, and the number of columns the
  // "width" of the matrix
  El::Int height = static_cast<El::Int>(numColumns());
  El::Int width  = static_cast<El::Int>(numRows());
  El::Int root   = static_cast<El::Int>(master_rank);
  El::DistMatrix<double> scratch(height, width, grid, root);

  // Columns of the matrix should be distributed round-robin; each
  // element of a column should be on the same process. The
  // distribution should satisfy two invariants:
  //
  // (1) Each "global row" should have the same rank in its column
  // communicator
  //
  // (2) For the scratch matrix global column index j, the process (on
  // the row communicator) that owns j should be process rank (j %
  // d_num_procs).
  //
  CAROM_ASSERT(scratch.RowOwner(0) == scratch.RowOwner(1));
  int my_rank;
  CAROM_ASSERT(MPI_Comm_rank(comm, &my_rank) == MPI_SUCCESS);
  El::Int rank_as_row = static_cast<El::Int>(my_rank);
  CAROM_ASSERT(scratch.ColOwner(rank_as_row) == rank_as_row);

  // Set up work matrices
  El::DistMatrix<double> householder_scalars(grid);
  El::DistMatrix<double> diagonal(grid);
  El::DistPermutation perm(grid);
  El::QRCtrl<double> ctrl;

  // Build the copy of the matrix element by element, mapping local
  // indices to global indices. The Elemental::DistMatrix indices are
  // the transpose of the Matrix indices. The El::Grid object should
  // be constructed so that each column of the scratch matrix is owned
  // by a single process.
  //
  // If the matrix is unbalanced, matrix elements need to be
  // redistributed among processes; then, for the purposes of
  // computing the QR decomposition only, we redistribute matrix
  // elements in the scratch matrix. First, we compute global to
  // (process rank, local) index map. The mapping chosen is not
  // terribly performant -- it may do more data movement than
  // necessary when assigning matrix elements -- but is easy to
  // implement.

  int *row_offset = new int[d_num_procs + 1];
  row_offset[d_num_procs + 1] = num_total_rows;

  row_offset[my_rank] = d_num_rows;
  const int send_recv_count = 1;
  CAROM_ASSERT(MPI_Allgather(MPI_IN_PLACE,
			     send_recv_count,
			     MPI_INT,
			     row_offset,
			     send_recv_count,
			     MPI_INT,
			     comm) == MPI_SUCCESS);

  for (size_t i = d_num_procs - 1; i >= 0; i--) {
    row_offset[i] = row_offset[i + 1] - row_offset[i];
  }
  CAROM_ASSERT(row_offset[0] == 0);

  // Use the computed (local <--> global) index map to assign elements
  // to the scratch matrix
  for (int row = row_offset[my_rank]; row < row_offset[my_rank + 1]; row++) {
    El::Int el_loc_col = static_cast<El::Int>(row - row_offset[my_rank]);
    El::Int el_global_col = scratch.GlobalCol(el_loc_col);
    for (int col = 0; col < d_num_cols; col++) {
      El::Int el_loc_row = static_cast<El::Int>(col);
      El::Int el_global_row = scratch.GlobalRow(el_loc_row);
      scratch.Set(el_global_row, el_global_col, this->item(row, col));
    }
  }

  // After transferring the data over to Elemental's native data
  // types, we compute the QRCP.
  El::QR(scratch, householder_scalars, diagonal, perm); // add ctrl if needed

  // Then, we transfer the pivots into the pivot array
  // stored redundantly on each process.
  for (size_t i = 0; i < pivots_requested; i++) {
    El::Int el_i = static_cast<El::Int>(i);
    El::Int el_perm_i = perm.Image(el_i);

    // The permutation is computed in terms of global row indices of
    // the Matrix, so we need to compute the rank that owns the column
    // of scratch (remember: scratch is the transpose of this Matrix)
    // corresponding to this global index
    int global_row_index = static_cast<int>(el_perm_i);
    int rank;
    for (rank = 0; rank < d_num_procs; rank++) {
      bool is_at_or_above_lower_bound = (global_row_index >= row_offset[rank]);
      bool is_below_upper_bound = (global_row_index < row_offset[rank + 1]);
      if (is_at_or_above_lower_bound && is_below_upper_bound) {
	row_pivot_owner[i] = rank;
	break;
      }
    }

    // The local row index is then the global row index minus the
    // row_offset computed above in the simple (local <--> global) index
    // map computed above
    row_pivot[i] = global_row_index - row_offset[rank];
  }

  // Free arrays
  delete [] row_offset;
#else
  CAROM_ASSERT(false);
#endif
}


} // end namespace CAROM
