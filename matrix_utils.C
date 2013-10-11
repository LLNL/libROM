#include <assert.h>
#include <mpi.h>

// Multiplies 2 matrices which are completely local to the calling process.
// Returns a matrix which is completely local to that process.
double*
LocalMatLocalMatMult(
   const double* a,
   int num_a_rows,
   int num_a_cols,
   const double* b,
   int num_b_rows,
   int num_b_cols)
{
   assert(num_a_cols == num_b_rows);
   double *result = new double [num_a_rows*num_b_cols];
   int result_idx = 0;
   for (int a_row = 0; a_row < num_a_rows; ++a_row) {
      for (int b_col = 0; b_col < num_b_cols; ++b_col) {
         int a_idx = a_row*num_a_cols;
         int b_idx = b_col;
         result[result_idx] = 0.0;
         for (int entry = 0; entry < num_a_cols; ++entry) {
            result[result_idx] += a[a_idx]*b[b_idx];
            ++a_idx;
            b_idx += num_b_cols;
         }
         ++result_idx;
      }
   }
   return result;
}

// Multiplies the transpose of a matrix, a, with another matrix, b, both of
// which are completely local to the calling process.  Returns a matrix which
// is completely local to that process.
double*
LocalMatTransposeLocalMatMult(
   const double* a,
   int num_a_rows,
   int num_a_cols,
   const double* b,
   int num_b_rows,
   int num_b_cols)
{
   assert(num_a_rows == num_b_rows);
   double *result = new double [num_a_cols*num_b_cols];
   int result_idx = 0;
   for (int a_col = 0; a_col < num_a_cols; ++a_col) {
      for (int b_col = 0; b_col < num_b_cols; ++b_col) {
         int a_idx = a_col;
         int b_idx = b_col;
         result[result_idx] = 0.0;
         for (int entry = 0; entry < num_a_rows; ++entry) {
            result[result_idx] += a[a_idx]*b[b_idx];
            a_idx += num_a_cols;
            b_idx += num_b_cols;
         }
         ++result_idx;
      }
   }
   return result;
}

// Multiplies a matrix whose rows are distributed across multiple processes, a,
// with a matrix which is completely local to the calling process, b.  Returns
// the part of the distributed result local to the calling process.
double*
DistributedMatLocalMatMult(
   const double* local_a,
   int num_local_a_rows,
   int num_a_cols,
   const double* b,
   int num_b_rows,
   int num_b_cols)
{
   return LocalMatLocalMatMult(local_a,
                               num_local_a_rows,
                               num_a_cols,
                               b,
                               num_b_rows,
                               num_b_cols);
}

// Multiplies the transpose of a matrix whose rows are distributed across
// multiple processes, a, with another matrix whose rows are distributed across
// multiple processes, b.  Returns the part of the distributed result local to
// the calling process.
double*
DistributedMatTransposeDistributedMatMult(
   const double* local_a,
   int num_local_a_rows,
   int num_a_cols,
   const double* local_b,
   int num_local_b_rows,
   int num_b_cols,
   int num_procs)
{
   assert(num_local_a_rows == num_local_b_rows);
   int new_mat_size = num_a_cols*num_b_cols;
   double *local_result = new double [new_mat_size];
   int result_idx = 0;
   for (int a_col = 0; a_col < num_a_cols; ++a_col) {
      for (int b_col = 0; b_col < num_b_cols; ++b_col) {
         int a_idx = a_col;
         int b_idx = b_col;
         local_result[result_idx] = 0.0;
         for (int entry = 0; entry < num_local_a_rows; ++entry) {
            local_result[result_idx] += local_a[a_idx]*local_b[b_idx];
            a_idx += num_a_cols;
            b_idx += num_b_cols;
         }
         ++result_idx;
      }
   }
   double* result;
   if (num_procs > 1) {
      result = new double [new_mat_size];
      MPI_Allreduce(local_result, result, new_mat_size, MPI_DOUBLE,
                    MPI_SUM, MPI_COMM_WORLD);
      delete local_result;
   }
   else {
      result = local_result;
   }
   return result;
}

// Multiplies the transpose of a matrix which is completely local to the
// calling process with another matrix whose rows are distributed across
// multiple processes.  Returns the part of the distributed result local to the
// calling process.
double*
LocalMatTransposeDistributedMatMult(
   const double* a,
   int num_a_rows,
   int num_a_cols,
   const double* local_b,
   int num_local_b_rows,
   int num_local_b_cols,
   int num_procs,
   int rank)
{
   assert(num_a_rows == num_local_b_rows);
   int new_mat_size = num_a_cols*num_local_b_cols;
   double *local_result = new double [new_mat_size];
   int result_idx = 0;
   for (int a_col = 0; a_col < num_a_cols; ++a_col) {
      for (int b_col = 0; b_col < num_local_b_cols; ++b_col) {
         int a_idx = a_col + num_a_rows*num_a_cols*rank;
         int b_idx = b_col;
         local_result[result_idx] = 0.0;
         for (int entry = 0; entry < num_a_rows; ++entry) {
            local_result[result_idx] += a[a_idx]*local_b[b_idx];
            a_idx += num_a_cols;
            b_idx += num_local_b_cols;
         }
         ++result_idx;
      }
   }
   double* result;
   if (num_procs > 1) {
      result = new double [new_mat_size];
      MPI_Allreduce(local_result, result, new_mat_size, MPI_DOUBLE,
                    MPI_SUM, MPI_COMM_WORLD);
      delete local_result;
   }
   else {
      result = local_result;
   }
   return result;
}
