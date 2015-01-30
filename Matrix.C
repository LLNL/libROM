/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2015 Lawrence Livermore National Security, LLC
 * Description: A simple, parallel Matrix class with the utility needed to
 *              support the basis generation methods of this library.  A
 *              distributed Matrix has its rows distributed across processors.
 *
 *****************************************************************************/

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
         result->item(this_row, other_col) = 0.0;
         for (int entry = 0; entry < d_num_cols; ++entry) {
            result->item(this_row, other_col) +=
               item(this_row, entry)*other.item(entry, other_col);
         }
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
      result->item(this_row) = 0.0;
      for (int entry = 0; entry < d_num_cols; ++entry) {
         result->item(this_row) += item(this_row, entry)*other.item(entry);
      }
   }
   return result;
}

Matrix*
Matrix::transposeMult(
   const Matrix& other) const
{
   CAROM_ASSERT(distributed() == other.distributed());
   CAROM_ASSERT(numRows() == other.numRows());
   if (!d_distributed) {
      Matrix* result = new Matrix(d_num_cols,
                                  other.d_num_cols,
                                  false);
      for (int this_col = 0; this_col < d_num_cols; ++this_col) {
         for (int other_col = 0; other_col < other.d_num_cols; ++other_col) {
            result->item(this_col, other_col) = 0.0;
            for (int entry = 0; entry < d_num_rows; ++entry) {
               result->item(this_col, other_col) +=
                  item(entry, this_col)*other.item(entry, other_col);
            }
         }
      }
      return result;
   }
   else {
      int new_mat_size = d_num_cols*other.d_num_cols;
      Matrix* local_result = new Matrix(d_num_cols,
                                        other.d_num_cols,
                                        false);
      for (int this_col = 0; this_col < d_num_cols; ++this_col) {
         for (int other_col = 0; other_col < other.d_num_cols; ++other_col) {
            local_result->item(this_col, other_col) = 0.0;
            for (int entry = 0; entry < d_num_rows; ++entry) {
               local_result->item(this_col, other_col) +=
                  item(entry, this_col)*other.item(entry, other_col);
            }
         }
      }
      Matrix* result;
      int mpi_init;
      MPI_Initialized(&mpi_init);
      int num_procs;
      if (mpi_init) {
         MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
      }
      else {
         num_procs = 1;
      }
      if (num_procs > 1) {
         result = new Matrix(d_num_cols,
                             other.d_num_cols,
                             false);
         MPI_Allreduce(local_result->d_mat, result->d_mat, new_mat_size,
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         delete local_result;
      }
      else {
         result = local_result;
      }
      return result;
   }
}

Vector*
Matrix::transposeMult(
   const Vector& other) const
{
   CAROM_ASSERT(distributed() && other.distributed());
   CAROM_ASSERT(numRows() == other.dim());
   Vector* local_result = new Vector(d_num_cols, false);
   int result_dim = 0;
   for (int this_col = 0; this_col < d_num_cols; ++this_col) {
      int other_dim = 0;
      local_result->item(result_dim) = 0.0;
      for (int entry = 0; entry < d_num_rows; ++entry) {
         local_result->item(result_dim) +=
            item(entry, this_col)*other.item(other_dim);
         ++other_dim;
      }
      ++result_dim;
   }
   Vector* result;
   int mpi_init;
   MPI_Initialized(&mpi_init);
   int num_procs;
   if (mpi_init) {
      MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   }
   else {
      num_procs = 1;
   }
   if (num_procs > 1) {
      result = new Vector(d_num_cols, false);
      MPI_Allreduce(&local_result->item(0), &result->item(0), d_num_cols,
                    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      delete local_result;
   }
   else {
      result = local_result;
   }
   return result;
}

}
