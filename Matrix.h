/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: A simple, parallel Matrix class with the utility needed to
 *              support the basis generation methods of this library.  A
 *              distributed Matrix has its rows distributed across processors.
 *
 *****************************************************************************/

#ifndef included_Matrix_h
#define included_Matrix_h

#include "Vector.h"

namespace CAROM {

/**
 * A simple matrix class in which the rows may be distributed across multuple
 * processes.  This class supports only the basic operations that are needed by
 * the SVD library.
 */
class Matrix
{
   public:
      /** Constructor creating a Matrix with uninitialized values.
       *
       * @pre num_rows > 0
       * @pre num_cols > 0
       *
       * @param[in] num_rows When undistributed, the total number of rows of
       *                     the Matrix.  When distributed, the part of the
       *                     total number of rows of the Matrix on this
       *                     processor.
       * @param[in] num_cols The total number of columns of the Matrix.
       * @param[in] distributed If true the rows of the Matrix are spread over
       *                        all processors.
       */
      Matrix(
         int num_rows,
         int num_cols,
         bool distributed);

      /** Constructor creating a Matrix with uninitialized values.
       *
       * @pre mat != 0
       * @pre num_rows > 0
       * @pre num_cols > 0
       *
       * @param[in] mat The initial values of the Matrix.
       * @param[in] num_rows When undistributed, the total number of rows of
       *                     the Matrix.  When distributed, the part of the
       *                     total number of rows of the Matrix on this
       *                     processor.
       * @param[in] num_cols The total number of columns of the Matrix.
       * @param[in] distributed If true the rows of the Matrix are spread over
       *                        all processors.
       */
      Matrix(
         const double* mat,
         int num_rows,
         int num_cols,
         bool distributed);

      /**
       * @brief Copy constructor.
       *
       * @param[in] other The Matrix to copy.
       */
      Matrix(
         const Matrix& other);

      /**
       * @brief Destructor.
       */
      ~Matrix();

      /**
       * @brief Assignment operator.
       *
       * @param[in] rhs The Matrix to assign to this.
       *
       * @return This after rhs has been assigned to it.
       */
      Matrix&
      operator = (
         const Matrix& rhs);

      /**
       * @brief Returns true if the Matrix is distributed.
       *
       * @return True if the Matrix is distributed.
       */
      bool
      distributed() const
      {
         return d_distributed;
      }

      /**
       * @brief Returns the number of rows of the Matrix on this processor.
       *
       * @return The number of rows of the Matrix on this processor.
       */
      int
      numRows() const
      {
         return d_num_rows;
      }

      /**
       * @brief Returns the number of columns in the Matrix.
       *
       * This method will return the same value from each processor.
       *
       * @return The number of columns of the Matrix.
       */
      int
      numColumns() const
      {
         return d_num_cols;
      }

      /**
       * @brief Multiplies this Matrix with other and returns the product,
       * reference version.
       *
       * Supports multiplication of two undistributed matrices returning an
       * undistributed Matrix, and multiplication of a distributed Matrix with
       * an undistributed Matrix returning a distributed Matrix.
       *
       * @pre !other.d_distributed
       * @pre d_num_cols == other.d_num_rows
       *
       * @param[in] other The Matrix to multiply with this.
       *
       * @return The product Matrix.
       */
      Matrix*
      mult(
         const Matrix& other) const;

      /**
       * @brief Multiplies this Matrix with other and returns the product,
       * pointer version.
       *
       * Supports multiplication of two undistributed matrices returning an
       * undistributed Matrix, and multiplication of a distributed Matrix with
       * an undistributed Matrix returning a distributed Matrix.
       *
       * @pre other != 0
       * @pre !other.d_distributed
       * @pre d_num_cols == other.d_num_rows
       *
       * @param[in] other The Matrix to multiply with this.
       *
       * @return The product Matrix.
       */
      Matrix*
      mult(
         const Matrix* const other) const
      {
         CAROM_ASSERT(other != 0);
         return mult(*other);
      }

      /**
       * @brief Multiplies this Matrix with other and returns the product,
       * reference version.
       *
       * Supports multiplication of a distributed Matrix and an undistributed
       * Vector returning a distributed Vector.
       *
       * @pre d_distributed && !other.distributed()
       * @pre d_num_cols == other.dim()
       *
       * @param[in] other The Vector to multiply with this.
       *
       * @return The product Vector.
       */
      Vector*
      mult(
         const Vector& other) const;

      /**
       * @brief Multiplies this Matrix with other and returns the product,
       * pointer version.
       *
       * Supports multiplication of a distributed Matrix and an undistributed
       * Vector returning a distributed Vector.
       *
       * @pre other != 0
       * @pre d_distributed && !other.distributed()
       * @pre d_num_cols == other.dim()
       *
       * @param[in] other The Vector to multiply with this.
       *
       * @return The product Vector.
       */
      Vector*
      mult(
         const Vector* const other) const
      {
         CAROM_ASSERT(other != 0);
         return mult(*other);
      }

      /**
       * @brief Multiplies the transpose of this Matrix with other and returns
       * the product, reference version.
       *
       * Supports multiplication of two undistributed matrices returning an
       * undistributed Matrix or two distributed matrices returning an
       * undistributed Matrix.
       *
       * @pre d_distributed == other.d_distributed
       * @pre d_num_rows == other.d_num_rows
       *
       * @param[in] other The Matrix to multiply with this.
       *
       * @return The product Matrix.
       */
      Matrix*
      transposeMult(
         const Matrix& other) const;

      /**
       * @brief Multiplies the transpose of this Matrix with other and returns
       * the product, pointer version.
       *
       * Supports multiplication of two undistributed matrices returning an
       * undistributed Matrix or two distributed matrices returning an
       * undistributed Matrix.
       *
       * @pre other != 0
       * @pre d_distributed == other.d_distributed
       * @pre d_num_rows == other.d_num_rows
       *
       * @param[in] other The Matrix to multiply with this.
       *
       * @return The product Matrix.
       */
      Matrix*
      transposeMult(
         const Matrix* const other) const
      {
         CAROM_ASSERT(other != 0);
         return transposeMult(*other);
      }

      /**
       * @brief Multiplies the transpose of this Matrix with other and returns
       * the product, reference version.
       *
       * Supports multiplication of a distributed Matrix and a distributed
       * Vector returning an undistributed Vector.
       *
       * @pre d_distributed && other.distributed()
       * @pre d_num_rows == other.dim();
       *
       * @param[in] other The Vector to multiply with this.
       *
       * @return The product Vector.
       */
      Vector*
      transposeMult(
         const Vector& other) const;

      /**
       * @brief Multiplies the transpose of this Matrix with other and returns
       * the product, pointer version.
       *
       * Supports multiplication of a distributed Matrix and a distributed
       * Vector returning an undistributed Vector.
       *
       * @pre other != 0
       * @pre d_distributed && other.distributed()
       * @pre d_num_rows == other.dim();
       *
       * @param[in] other The Vector to multiply with this.
       *
       * @return The product Vector.
       */
      Vector*
      transposeMult(
         const Vector* const other) const
      {
         CAROM_ASSERT(other != 0);
         return transposeMult(*other);
      }

      /**
       * @brief Const Matrix member access.
       *
       * @pre (0 <= row) && (row < d_num_rows)
       * @pre (0 <= col) && (col < d_num_cols)
       *
       * @param[in] row The row of the Matrix value on this processor
       *                requested.
       * @param[in] col The column of the Matrix value requested.
       */
      const double&
      item(
         int row,
         int col) const
      {
         CAROM_ASSERT((0 <= row) && (row < d_num_rows));
         CAROM_ASSERT((0 <= col) && (col < d_num_cols));
         return d_mat[row*d_num_cols+col];
      }

      /**
       * @brief Non-const Matrix member access.
       *
       * Allows constructs of the form mat[i, j] = val;
       *
       * @pre (0 <= row) && (row < d_num_rows)
       * @pre (0 <= col) && (col < d_num_cols)
       *
       * @param[in] row The row of the Matrix value on this processor
       *                requested.
       * @param[in] col The column of the Matrix value requested.
       */
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
      /**
       * @brief Default constructor is not implemented.
       */
      Matrix();

      /**
       * @brief The storage for the Matrix's values on this processor.
       */
      double* d_mat;

      /**
       * @brief The rows in the Matrix that are on this processor.
       */
      int d_num_rows;

      /**
       * @brief The number of columns in the Matrix.
       *
       * For distributed matrices the number of columns is the same on all
       * processors.
       */
      int d_num_cols;

      /**
       * @brief If true, the Matrix's rows are distributed over all processors.
       *
       * Each processor does not need to hold the same number of rows.
       */
      bool d_distributed;
};

}

#endif
