#ifndef included_matrix_h
#define included_matrix_h

#include "vector.h"

namespace CAROM {

// A simple matrix class in which the rows may be distributed evenly across
// multiple processes.  Supports only the basic multiplication operations that
// we need in the CAROM project.
class Matrix
{
   public:
      // Constructor in which matrix manages storage.
      Matrix(
         int num_rows,
         int num_cols,
         bool distributed,
         int rank,
         int num_procs);

      // Copy constructor.
      Matrix(
         const Matrix& other);

      // Destructor.
      ~Matrix();

      // Returns true if matrix is distributed.
      bool
      distributed() const
      {
         return d_distributed;
      }

      // Returns the number of rows in the matrix.
      int
      numRows() const
      {
         return d_num_rows;
      }

      // Returns the number of columns in the matrix.
      int
      numColumns() const
      {
         return d_num_cols;
      }

      // Multiplies this matrix with other and returns the product.  Supports
      // multiplication of two undistributed matrices returning an
      // undistributed matrix, and multiplication of a distributed matrix with
      // an undistributed matrix returning a distributed matrix.
      Matrix*
      Mult(
         const Matrix& other) const;

      // Multiplies this matrix with other and returns the product.  Supports
      // multiplication of a distributed matrix and an undistributed vector
      // vector returning a distributed vector.
      Vector*
      Mult(
         const Vector& other) const;

      // Multiplies the transpose of this matrix with other and returns the
      // product.  Supports multiplication of two undistributed matrices
      // returning an undistributed matrix, two distributed matrices returning
      // an undistributed matrix, and an undistributed matrix with a
      // distributed matrix returning an undistributed matrix.
      Matrix*
      TransposeMult(
         const Matrix& other) const;

      // Multiplies the transpose of this matrix with other and returns the
      // product.  Supports multiplication of a distributed matrix and a
      // distributed vector returning an undistributed vector.
      Vector*
      TransposeMult(
         const Vector& other) const;

      // Const matrix member access.
      const double&
      item(
         int row,
         int col) const
      {
         CAROM_ASSERT((0 <= row) && (row < d_num_rows));
         CAROM_ASSERT((0 <= col) && (col < d_num_cols));
         return d_mat[row*d_num_cols+col];
      }

      // Non-const matrix member access.
      double&
      item(
         int row,
         int col)
      {
         CAROM_ASSERT((0 <= row) && (row < d_num_rows));
         CAROM_ASSERT((0 <= col) && (col < d_num_cols));
         return d_mat[row*d_num_cols+col];
      }

   private:
      // Default constructor is not implemented.
      Matrix();

      // Assignment operator is not implemented.
      Matrix&
      operator = (
         const Matrix& rhs);

      // The storage for the matrix.
      double* d_mat;

      // The number of rows in the matrix.
      int d_num_rows;

      // The number of columns in the matrix.
      int d_num_cols;

      // If true, the matrix's rows are distributed over all processors.  Each
      // processor holds the same number of rows.
      bool d_distributed;

      // The MPI rank of the process owning this object.
      int d_rank;

      // The number of MPI processes.
      int d_num_procs;
};

}

#endif
