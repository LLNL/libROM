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

      // Constructor in which matrix is given storage managed by some other
      // entity.
      Matrix(
         double* mat,
         int num_rows,
         int num_cols,
         bool distributed,
         int rank,
         int num_procs);

      // Destructor.
      ~Matrix();

      // Multiplies this matrix with other and returns the product.  Supports
      // multiplication of two undistributed matrices returning an
      // undistributed matrix, and multiplication of a distributed matrix with
      // an undistributed matrix returning a distributed matrix.
      Matrix*
      Mult(const Matrix& other) const;

      // Multiplies the transpose of this matrix with other and returns the
      // product.  Supports multiplication of two undistributed matrices
      // returning an undistributed matrix, two distributed matrices returning
      // an undistributed matrix, and an undistributed matrix with a
      // distributed matrix returning an undistributed matrix.
      Matrix*
      TransposeMult(const Matrix& other) const;

      // Multiplies the transpose of this matrix with other and returns the
      // product.  Supports multiplication of a distributed matrix and a
      // distributed vector returning an undistributed vector.
      Vector*
      TransposeMult(const Vector& other) const;

      // Const matrix member access.
      const double&
      item(const int row, const int col) const
      {
#ifdef DEBUG
         assert((0 <= row) && (row < d_num_rows));
         assert((0 <= col) && (col < d_num_cols));
#endif
         return d_mat[row*d_num_cols+col];
      }

      // Non-const matrix member access.
      double&
      item(const int row, const int col)
      {
#ifdef DEBUG
         assert((0 <= row) && (row < d_num_rows));
         assert((0 <= col) && (col < d_num_cols));
#endif
         return d_mat[row*d_num_cols+col];
      }

   private:
      friend class incremental_svd;
      friend class static_svd;

      // Default constructor is not implemented.
      Matrix();

      // Copy constructor is not implemented.
      Matrix(
         const Matrix& other);

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
      int d_distributed;

      // The MPI rank of the process owning this object.
      int d_rank;

      // The number of MPI processes.
      int d_num_procs;

      // True if the object manages its storage.
      bool d_manages_storage;
};

}

#endif
