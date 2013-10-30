#include "matrix.h"
#include <mpi.h>

namespace CAROM {

Matrix::Matrix(
   int num_rows,
   int num_cols,
   bool distributed,
   int rank,
   int num_procs) :
   d_num_rows(num_rows),
   d_num_cols(num_cols),
   d_distributed(distributed),
   d_rank(rank),
   d_num_procs(num_procs),
   d_manages_storage(true)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(num_rows > 0);
   assert(num_cols > 0);
   assert(rank < num_procs);
#endif
   d_mat = new double [num_cols*num_rows];
}

Matrix::Matrix(
   double* mat,
   int num_rows,
   int num_cols,
   bool distributed,
   int rank,
   int num_procs) :
   d_num_rows(num_rows),
   d_num_cols(num_cols),
   d_distributed(distributed),
   d_rank(rank),
   d_num_procs(num_procs),
   d_manages_storage(false)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(num_rows > 0);
   assert(num_cols > 0);
   assert(rank < num_procs);
#endif
   d_mat = mat;
}

Matrix::~Matrix()
{
   if (d_manages_storage) {
      delete [] d_mat;
   }
}

Matrix*
Matrix::Mult(
   const Matrix& other) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!other.d_distributed);
   assert(d_num_cols == other.d_num_rows);
#endif
   Matrix* result = new Matrix(d_num_rows,
                               other.d_num_cols,
                               d_distributed,
                               d_rank,
                               d_num_procs);
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

Matrix*
Matrix::TransposeMult(
   const Matrix& other) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert((d_distributed == other.d_distributed) ||
          (!d_distributed && other.d_distributed));
   assert((!d_distributed && other.d_distributed) ||
          (d_num_rows == other.d_num_rows));
#endif
   if (!d_distributed && !other.d_distributed) {
      Matrix* result = new Matrix(d_num_cols,
                                  other.d_num_cols,
                                  false,
                                  d_rank,
                                  d_num_procs);
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
   else if (d_distributed && other.d_distributed) {
      int new_mat_size = d_num_cols*other.d_num_cols;
      Matrix* local_result = new Matrix(d_num_cols,
                                        other.d_num_cols,
                                        false,
                                        d_rank,
                                        d_num_procs);
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
      if (d_num_procs > 1) {
         result = new Matrix(d_num_cols,
                             other.d_num_cols,
                             false,
                             d_rank,
                             d_num_procs);
         MPI_Allreduce(local_result->d_mat, result->d_mat, new_mat_size,
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         delete local_result;
      }
      else {
         result = local_result;
      }
      return result;
   }
   else {
      int new_mat_size = d_num_cols*other.d_num_cols;
      Matrix* local_result = new Matrix(d_num_cols,
                                        other.d_num_cols,
                                        false,
                                        d_rank,
                                        d_num_procs);
      for (int this_col = 0; this_col < d_num_cols; ++this_col) {
         for (int other_col = 0; other_col < other.d_num_cols; ++other_col) {
            int this_row = other.d_num_rows*d_rank;
            local_result->item(this_col, other_col) = 0.0;
            for (int entry = 0; entry < other.d_num_rows; ++entry) {
               local_result->item(this_col, other_col) +=
                  item(this_row, this_col)*other.item(entry, other_col);
               ++this_row;
            }
         }
      }
      Matrix* result;
      if (d_num_procs > 1) {
         result = new Matrix(d_num_cols,
                             other.d_num_cols,
                             false,
                             d_rank,
                             d_num_procs);
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
Matrix::TransposeMult(const Vector& other) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(d_distributed && other.d_distributed);
   assert(d_num_rows == other.d_dim);
#endif
   Vector* local_result = new Vector(d_num_cols, false, d_rank, d_num_procs);
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
   if (d_num_procs > 1) {
      result = new Vector(d_num_cols, false, d_rank, d_num_procs);
      MPI_Allreduce(local_result->d_vec, result->d_vec, d_num_cols,
                    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      delete local_result;
   }
   else {
      result = local_result;
   }
   return result;
}

}
