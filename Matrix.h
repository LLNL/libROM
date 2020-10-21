/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: A simple, parallel Matrix class with the utility needed to
//              support the basis generation methods of this library.  A
//              distributed Matrix has its rows distributed across processors.

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
      /** Empty Constructor */
      Matrix();

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
       * @param[in] copy_data If true the matrix allocates is own storage and
       *                      copies the contents of mat into its own storage.
       *                      Otherwise it uses mat as its storage.
       */
      Matrix(
         double* mat,
         int num_rows,
         int num_cols,
         bool distributed,
         bool copy_data = true);

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
       * @brief Assignment operator.
       *
       * @param[in] a constant value
       *
       * @return This after filling all the data with a constant value
       */
      Matrix&
      operator = (
         const double a);

      /**
       * @brief Addition operator.
       *
       * @param[in] rhs The Matrix to add to this.
       *
       * @return This after rhs has been added to it.
       */
      Matrix&
      operator += (
         const Matrix& rhs);

      /**
       * @brief Subtraction operator.
       *
       * @param[in] rhs The Matrix to subtract to this.
       *
       * @return This after rhs has been subtracted to it.
       */
      Matrix&
      operator -= (
         const Matrix& rhs);

      /**
       * @brief Sets the number of rows and columns of the matrix and
       * reallocates storage if needed.
       *
       * @param[in] num_rows New number of rows
       * @param[in] num_cols New number of cols
       */
      void
      setSize(
         int num_rows,
         int num_cols)
      {
         int new_size = num_rows*num_cols;
         if (new_size > d_alloc_size) {
            if (!d_owns_data) {
               CAROM_ERROR("Can not reallocate externally owned storage.");
            }
            if (d_mat) {
               delete [] d_mat;
            }
            d_mat = new double [new_size];
            d_alloc_size = new_size;
         }
         d_num_rows = num_rows;
         d_num_cols = num_cols;
      }

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
       * @brief Returns true if rows of matrix are load-balanced.
       *
       */
      bool balanced() const;

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
       * @pre !other.distributed()
       * @pre numColumns() == other.numRows()
       *
       * @param[in] other The Matrix to multiply with this.
       *
       * @return The product Matrix.
       */
      Matrix*
      mult(
         const Matrix& other) const
      {
        Matrix* result = 0;
        mult(other, result);
        return result;
      }

      /**
       * @brief Multiplies this Matrix with other and returns the product,
       * pointer version.
       *
       * Supports multiplication of two undistributed matrices returning an
       * undistributed Matrix, and multiplication of a distributed Matrix with
       * an undistributed Matrix returning a distributed Matrix.
       *
       * @pre other != 0
       * @pre !other->distributed()
       * @pre numColumns() == other->numRows()
       *
       * @param[in] other The Matrix to multiply with this.
       *
       * @return The product Matrix.
       */
      Matrix*
      mult(
         const Matrix* other) const
      {
         CAROM_ASSERT(other != 0);
         return mult(*other);
      }

      /**
       * @brief Multiplies this Matrix with other and fills result with the
       * answer.
       *
       * Supports multiplication of two undistributed matrices resulting in an
       * undistributed Matrix, and multiplication of a distributed Matrix with
       * an undistributed Matrix resulting in a distributed Matrix.  If result
       * has not been allocated it will be, otherwise it will be sized
       * accordingly.
       *
       * @pre result == 0 || result->distributed() == distributed()
       * @pre !other.distributed()
       * @pre numColumns() == other.numRows()
       *
       * @param[in] other The Matrix to multiply with this.
       * @param[out] result The product Matrix.
       */
      void
      mult(
         const Matrix& other,
         Matrix*& result) const;

      /**
       * @brief Multiplies this Matrix with other and fills result with the
       * answer.
       *
       * Supports multiplication of two undistributed matrices resulting in an
       * undistributed Matrix, and multiplication of a distributed Matrix with
       * an undistributed Matrix resulting in a distributed Matrix.  Result
       * will be sized accordingly.
       *
       * @pre result.distributed() == distributed()
       * @pre !other.distributed()
       * @pre numColumns() == other.numRows()
       *
       * @param[in] other The Matrix to multiply with this.
       * @param[out] result The product Matrix.
       */
      void
      mult(
         const Matrix& other,
         Matrix& result) const;

      /**
       * @brief Multiplies this Matrix with other and returns the product,
       * reference version.
       *
       * Supports multiplication of an undistributed Matrix and Vector
       * returning an undistributed Vector, and multiplication of a distributed
       * Matrix and an undistributed Vector returning a distributed Vector.
       *
       * @pre !other.distributed()
       * @pre numColumns() == other.dim()
       *
       * @param[in] other The Vector to multiply with this.
       *
       * @return The product Vector.
       */
      Vector*
      mult(
         const Vector& other) const
      {
         Vector* result = 0;
         mult(other, result);
         return result;
      }

      /**
       * @brief Multiplies this Matrix with other and returns the product,
       * pointer version.
       *
       * Supports multiplication of an undistributed Matrix and Vector
       * returning an undistributed Vector, and multiplication of a distributed
       * Matrix and an undistributed Vector returning a distributed Vector.
       *
       * @pre other != 0
       * @pre !other->distributed()
       * @pre numColumns() == other->dim()
       *
       * @param[in] other The Vector to multiply with this.
       *
       * @return The product Vector.
       */
      Vector*
      mult(
         const Vector* other) const
      {
         CAROM_ASSERT(other != 0);
         return mult(*other);
      }

      /**
       * @brief Multiplies this Matrix with other and fills result with the
       * answer.
       *
       * Supports multiplication of an undistributed Matrix and Vector
       * resulting in an undistributed Vector, and multiplication of a
       * distributed Matrix and an undistributed Vector resulting in a
       * distributed Vector.  If result has not been allocated it will be,
       * otherwise it will be sized accordingly.
       *
       * @pre result == 0 || result->distributed() == distributed()
       * @pre !other.distributed()
       * @pre numColumns() == other.dim()
       *
       * @param[in] other The Vector to multiply with this.
       * @param[out] result The product Vector.
       */
      void
      mult(
         const Vector& other,
         Vector*& result) const;

      /**
       * @brief Multiplies this Matrix with other and fills result with the
       * answer.
       *
       * Supports multiplication of an undistributed Matrix and Vector
       * resulting in an undistributed Vector, and multiplication of a
       * distributed Matrix and an undistributed Vector resulting in a
       * distributed Vector.  Result will be sized accordingly.
       *
       * @pre result.distributed() == distributed()
       * @pre !other.distributed()
       * @pre numColumns() == other.dim()
       *
       * @param[in] other The Vector to multiply with this.
       * @param[out] result The product Vector.
       */
      void
      mult(
         const Vector& other,
         Vector& result) const;

        /**
        * @brief Multiplies a specified row of this Matrix with other
        * pointwise.
        *
        * Only supports multiplication of an undistributed Matrix and Vector
        * resulting in an undistributed Vector.
        * Result will be sized accordingly.
        *
        * @pre !result.distributed()
        * @pre !distributed()
        * @pre !other.distributed()
        * @pre numColumns() == other.dim()
        *
        * @param[in] other The Vector to multiply with this.
        * @param[out] result The product Vector.
        */
      void
      pointwise_mult(
         int this_row,
         const Vector& other,
         Vector& result) const;

        /**
        * @brief Multiplies a specified row of this Matrix with other
        * pointwise. This modifies other.
        *
        * Only supports multiplication of an undistributed Matrix and Vector
        * resulting in an undistributed Vector.
        * Result will be sized accordingly.
        *
        * @pre !result.distributed()
        * @pre !distributed()
        * @pre !other.distributed()
        * @pre numColumns() == other.dim()
        *
        * @param[in] other The Vector to multiply with this.
        * @param[out] other The product Vector.
        */
      void
      pointwise_mult(
         int this_row,
         Vector& other) const;


      /**
       * @brief Computes a += this*b*c.
       *
       * Supports accumulation of the multiplication of an undistributed
       * Matrix and Vector into an undistributed Vector, and accumulation of
       * the multiplication of a distributed Matrix and an undistributed
       * Vector into a distributed Vector.
       *
       * @pre a.distributed() == distributed()
       * @pre !b->distributed()
       * @pre numColumns() == b.dim()
       * @pre numRows() = a.dim()
       *
       * @param[in,out] a The Vector to accumulate this*b into.
       * @param[in] b The Vector multiplied by this.
       * @param[in] c Scalar multiplication factor.
       */
      void
      multPlus(
         Vector& a,
         const Vector& b,
         double c) const;

      /**
       * @brief Multiplies the transpose of this Matrix with other and returns
       * the product, reference version.
       *
       * Supports multiplication of two undistributed matrices returning an
       * undistributed Matrix or two distributed matrices returning an
       * undistributed Matrix.
       *
       * @pre distributed() == other.distributed()
       * @pre numRows() == other.numRows()
       *
       * @param[in] other The Matrix to multiply with this.
       *
       * @return The product Matrix.
       */
      Matrix*
      transposeMult(
         const Matrix& other) const
      {
         Matrix* result = 0;
         transposeMult(other, result);
         return result;
      }

      /**
       * @brief Multiplies the transpose of this Matrix with other and returns
       * the product, pointer version.
       *
       * Supports multiplication of two undistributed matrices returning an
       * undistributed Matrix or two distributed matrices returning an
       * undistributed Matrix.
       *
       * @pre other != 0
       * @pre distributed() == other->distributed()
       * @pre numRows() == other->numRows()
       *
       * @param[in] other The Matrix to multiply with this.
       *
       * @return The product Matrix.
       */
      Matrix*
      transposeMult(
         const Matrix* other) const
      {
         CAROM_ASSERT(other != 0);
         return transposeMult(*other);
      }

      /**
       * @brief Multiplies the transpose of this Matrix with other and fills
       * result with the answer.
       *
       * Supports multiplication of two undistributed matrices or two
       * distributed matrices resulting in an undistributed Matrix.  If result
       * has not been allocated it will be, otherwise it will be sized
       * accordingly.
       *
       * @pre result == 0 || !result->distributed()
       * @pre distributed() == other.distributed()
       * @pre numRows() == other.numRows()
       *
       * @param[in] other The Matrix to multiply with this.
       * @param[out] result The product Matrix.
       */
      void
      transposeMult(
         const Matrix& other,
         Matrix*& result) const;

      /**
       * @brief Multiplies the transpose of this Matrix with other and fills
       * result with the answer.
       *
       * Supports multiplication of two undistributed matrices or two
       * distributed matrices resulting in an undistributed Matrix.  Result
       * will be sized accordingly.
       *
       * @pre !result.distributed()
       * @pre distributed() == other.distributed()
       * @pre numRows() == other.numRows()
       *
       * @param[in] other The Matrix to multiply with this.
       * @param[out] result The product Matrix.
       */
      void
      transposeMult(
         const Matrix& other,
         Matrix& result) const;

      /**
       * @brief Multiplies the transpose of this Matrix with other and returns
       * the product, reference version.
       *
       * Supports multiplication of an undistributed Matrix and an
       * undistributed Vector or a distributed Matrix and a distributed Vector
       * returning an undistributed Vector.
       *
       * @pre distributed() == other.distributed()
       * @pre numRows() == other.dim();
       *
       * @param[in] other The Vector to multiply with this.
       *
       * @return The product Vector.
       */
      Vector*
      transposeMult(
         const Vector& other) const
      {
         Vector* result = 0;
         transposeMult(other, result);
         return result;
      }

      /**
       * @brief Multiplies the transpose of this Matrix with other and returns
       * the product, pointer version.
       *
       * Supports multiplication of an undistributed Matrix and an
       * undistributed Vector or a distributed Matrix and a distributed Vector
       * returning an undistributed Vector.
       *
       * @pre other != 0
       * @pre distributed() == other->distributed()
       * @pre numRows() == other->dim();
       *
       * @param[in] other The Vector to multiply with this.
       *
       * @return The product Vector.
       */
      Vector*
      transposeMult(
         const Vector* other) const
      {
         CAROM_ASSERT(other != 0);
         return transposeMult(*other);
      }

      /**
       * @brief Multiplies the transpose of this Matrix with other and fills
       * result with the answer.
       *
       * Supports multiplication of an undistributed Matrix and an
       * undistributed Vector or a distributed Matrix and a distributed Vector
       * resulting in an undistributed Vector.  If result has not been allocated
       * it will be, otherwise it will be sized accordingly.
       *
       * @pre result == 0 || !result->distributed()
       * @pre distributed() == other.distributed()
       * @pre numRows() == other.dim();
       *
       * @param[in] other The Vector to multiply with this.
       * @param[out] result The product Vector.
       */
      void
      transposeMult(
         const Vector& other,
         Vector*& result) const;

      /**
       * @brief Multiplies the transpose of this Matrix with other and fills
       * result with the answer.
       *
       * Supports multiplication of an undistributed Matrix and an
       * undistributed Vector or a distributed Matrix and a distributed Vector
       * resulting in an undistributed Vector.  Result will be sized
       * accordingly.
       *
       * @pre !result.distributed()
       * @pre distributed() == other.distributed()
       * @pre numRows() == other.dim();
       *
       * @param[in] other The Vector to multiply with this.
       * @param[out] result The product Vector.
       */
      void
      transposeMult(
         const Vector& other,
         Vector& result) const;

      /**
       * @brief Computes and returns the inverse of this.
       *
       * @pre !distributed()
       * @pre numRows() == numColumns()
       *
       * @return The inverse of this.
       */
      Matrix*
      inverse() const
      {
         Matrix* result = 0;
         inverse(result);
         return result;
      }

      /**
       * @brief Computes and returns the inverse of this.
       *
       * If result has not been allocated it will be, otherwise it will be
       * sized accordingly.
       *
       * @pre result == 0 || (!result->distributed() &&
       *                      result->numRows() == numRows() &&
       *                      result->numColumns() == numColumns())
       * @pre !distributed()
       * @pre numRows() == numColumns()
       *
       * @param[out] result The inverse of this.
       */
      void
      inverse(
         Matrix*& result) const;

      /**
       * @brief Computes and returns the inverse of this.
       *
       * Result will be sized accordingly.
       *
       * @pre !result.distributed() && result.numRows() == numRows() &&
       *      result.numColumns() == numColumns()
       * @pre !distributed()
       * @pre numRows() == numColumns()
       *
       * @param[out] result The inverse of this.
       */
      void
      inverse(
         Matrix& result) const;

      /**
       * @brief Computes the inverse of this and stores result in this.
       *
       * @pre !distributed()
       * @pre numRows() == numColumns()
       */
      void
      inverse();

      /**
       * @brief Replaces this Matrix with its transpose (in place),
       * in the serial square case only.
       *
       * @pre !distributed()
       * @pre numRows() == numColumns()
       */
      void transpose();

      /**
       * @brief Compute the leading numColumns() column pivots from a
       * QR decomposition with column pivots (QRCP) of the transpose
       * of this.
       *
       * @pre !distributed()
       *
       * @param[out] row_pivot Array of leading column pivots
       * from QRCP of transpose of this Matrix, has length pivots_requested
       * @param[out] row_pivot_owner Array of process rank that owns
       * each pivot on the communicator owned by this Matrix.
       * @param[in]  number of pivots requested, must be less than or equal
       * to the number of rows of this Matrix.
       */
      void
      qrcp_pivots_transpose(int* row_pivot,
			    int* row_pivot_owner,
			    int  pivots_requested) const;

      /**
       * @brief Const Matrix member access. Matrix data is stored in
       * row-major format.
       *
       * @pre (0 <= row) && (row < numRows())
       * @pre (0 <= col) && (col < numColumns())
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
         CAROM_ASSERT((0 <= row) && (row < numRows()));
         CAROM_ASSERT((0 <= col) && (col < numColumns()));
         return d_mat[row*d_num_cols+col];
      }

      /**
       * @brief Non-const Matrix member access. Matrix data is stored
       * in row-major format.
       *
       * Allows constructs of the form mat[i, j] = val;
       *
       * @pre (0 <= row) && (row < numRows())
       * @pre (0 <= col) && (col < numColumns())
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
         CAROM_ASSERT((0 <= row) && (row < numRows()));
         CAROM_ASSERT((0 <= col) && (col < numColumns()));
         return d_mat[row*d_num_cols+col];
      }

      /**
       * @brief Const Matrix member access.
       *
       * @pre (0 <= row) && (row < numRows())
       * @pre (0 <= col) && (col < numColumns())
       *
       * @param[in] row The row of the Matrix value on this processor
       *                requested.
       * @param[in] col The column of the Matrix value requested.
       */
      const double& operator() (int row, int col) const
      {
         return item(row, col);
      }

      /**
       * @brief Non-const Matrix member access.
       *
       * Allows constructs of the form mat[i, j] = val;
       *
       * @pre (0 <= row) && (row < numRows())
       * @pre (0 <= col) && (col < numColumns())
       *
       * @param[in] row The row of the Matrix value on this processor
       *                requested.
       * @param[in] col The column of the Matrix value requested.
       */
      double& operator() (int row, int col)
      {
         return item(row, col);
      }

     /**
      * @brief print Matrix into (a) ascii file(s).
      *
      * @param[in] prefix The name of the prefix of the file name.
      *
      */
     void print(const char * prefix);

     /**
      * @brief write Matrix into (a) HDF file(s).
      *
      * @param[in] prefix The name of the prefix of the file name.
      *
      */
     void write(const std::string& base_file_name);

     /**
      * @brief read Matrix into (a) HDF file(s).
      *
      * @param[in] prefix The name of the prefix of the file name.
      *
      */
     void read(const std::string& base_file_name);

     double *getData() const
     {
       return d_mat;
     }
     
   private:
      /**
       * @brief Default constructor is not implemented.
       */
      // Matrix();

      /**
       * @brief Compute the leading numColumns() column pivots from a
       * QR decomposition with column pivots (QRCP) of the transpose
       * of this.
       *
       * @pre !distributed()
       *
       * @param[out] row_pivot Array of leading column pivots
       * from QRCP of transpose of this Matrix, has length pivots_requested
       * @param[out] row_pivot_owner Array of process rank that owns
       * each pivot on the communicator owned by this Matrix.
       * @param[in]  number of pivots requested, must be less than or equal
       * to the number of rows of this Matrix.
       */
      void
      qrcp_pivots_transpose_serial(int* row_pivot,
				   int* row_pivot_owner,
				   int  pivots_requested) const;

      /**
       * @brief Compute the leading column pivots from a QR
       * decomposition with column pivots (QRCP) of the transpose of
       * this Matrix, if it is distributed and balanced.
       *
       * @pre distributed()
       *
       * @param[out] row_pivot Array of leading column pivots
       * from QRCP of transpose of this Matrix, has length pivots_requested
       * @param[out] row_pivot_owner Array of process rank that owns
       * each pivot on the communicator owned by this Matrix.
       * @param[in]  number of pivots requested, must be less than or equal
       * to the number of rows of this Matrix.
       */
      void
      qrcp_pivots_transpose_distributed(int* row_pivot,
					int* row_pivot_owner,
					int  pivots_requested) const;

      void
      qrcp_pivots_transpose_distributed_scalapack(int* row_pivot,
					int* row_pivot_owner,
					int  pivots_requested) const;

      /**
       * @brief Compute the leading column pivots from a QR
       * decomposition with column pivots (QRCP) of the transpose of
       * this Matrix, if it is distributed and balanced. Prototype
       * using Elemental (requires C++11).
       *
       * @pre distributed()
       *
       * @param[out] row_pivot Array of leading column pivots
       * from QRCP of transpose of this Matrix, has length pivots_requested
       * @param[out] row_pivot_owner Array of process rank that owns
       * each pivot on the communicator owned by this Matrix.
       * @param[in]  number of pivots requested, must be less than or equal
       * to the number of rows of this Matrix.
       */
      void
      qrcp_pivots_transpose_distributed_elemental(int* row_pivot,
						  int* row_pivot_owner,
						  int  pivots_requested) const;

      /**
       * @brief Compute the leading column pivots from a QR
       * decomposition with column pivots (QRCP) of the transpose of
       * this Matrix, if it is distributed and balanced. Prototype
       * using Elemental (requires C++11).
       *
       * @pre distributed() && balanced()
       *
       * @param[out] row_pivot Array of leading column pivots
       * from QRCP of transpose of this Matrix, has length pivots_requested
       * @param[out] row_pivot_owner Array of process rank that owns
       * each pivot on the communicator owned by this Matrix.
       * @param[in]  number of pivots requested, must be less than or equal
       * to the number of rows of this Matrix.
       */
      void
      qrcp_pivots_transpose_distributed_elemental_balanced
	(int* row_pivot, int* row_pivot_owner, int pivots_requested) const;

      /**
       * @brief Compute the leading column pivots from a QR
       * decomposition with column pivots (QRCP) of the transpose of
       * this Matrix, if it is distributed and balanced. Prototype
       * using Elemental (requires C++11).
       *
       * @pre distributed() && !balanced()
       *
       * @param[out] row_pivot Array of leading column pivots
       * from QRCP of transpose of this Matrix, has length pivots_requested
       * @param[out] row_pivot_owner Array of process rank that owns
       * each pivot on the communicator owned by this Matrix.
       * @param[in]  number of pivots requested, must be less than or equal
       * to the number of rows of this Matrix.
       */
      void
      qrcp_pivots_transpose_distributed_elemental_unbalanced
	(int* row_pivot, int* row_pivot_owner, int pivots_requested) const;

      /**
       * @brief Computes the transposePseudoinverse of this.
       *
       * @pre !distributed()
       * @pre numRows() >= numColumns()
       *
       * Assumes this is full column rank; may fail if this is not
       * full column rank.
       */
      void transposePseudoinverse();
      friend void GNAT(const Matrix*,
		       const int,
		       int*,
		       int*,
		       Matrix&,
		       const int,
		       const int,
		       const int);

      friend void QDEIM(const Matrix*,
			const int,
			int*,
			int*,
			Matrix&,
			const int,
			const int,
			const int);

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
       * @brief The currently allocated size.
       *
       * d_num_row*d_num_cols <= d_alloc_size
       */
      int d_alloc_size;

      /**
       * @brief If true, the Matrix's rows are distributed over all processors.
       *
       * Each processor does not need to hold the same number of rows.
       */
      bool d_distributed;

      /**
       * @brief The number of processors being run on.
       */
      int d_num_procs;

      /**
       * @brief If true, this object owns its underlying data, d_mat, and
       * is responsible for its deletion.
       *
       * If d_owns_data is false, then the object may not reallocate d_mat.
       */
      bool d_owns_data;
};

    /**
     * @brief Computes the outer product of two Vectors, v and w.
     *
     * @post The number of rows in the returned matrix equals the dimension of v
     * @post The number of cols in the returned matrix equals the dimension of w,
     *       if w is not distributed. If w is distributed, the number of columns
     *       equals the total number of entries of w, summed over all processes.
     * @post If v is distributed, so is the returned matrix.
     *
     * @param[in] v The first vector in the outer product
     * @param[in] w The second vector in the outer product
     *
     * @return The product of v and the transpose of w.
     */
    // NOTE(goxberry@gmail.com; oxberry1@llnl.gov): Passing by value is
    // supposed to be more idiomatic C++11; consider passing by value
    // instead of passing by reference.
    Matrix outerProduct(const Vector &v, const Vector &w);

    /**
     * @brief Factory function to make a diagonal matrix with nonzero
     * entries as in its Vector argument. The rows of this diagonal
     * matrix are distributed in the same fashion as its Vector
     * argument.
     *
     * @param[in] v Vector of elements that will be on the diagonal of the
     *              output matrix.
     * @param[out] diagonalMatrix The diagonal matrix created by this function.
     *
     * @post diagonalMatrix.distributed() == v.distributed()
     * @post diagonalMatrix.numRows() == v.dim()
     */
    Matrix DiagonalMatrixFactory(const Vector &v);

    /**
     * @brief Factory function to make an identity matrix with rows
     * distributed in the same fashion as its Vector argument. This function
     * is provided for convenience due to the ubiquity of identity matrices
     * (and operators) in mathematics.
     *
     * @param[in] v Vector of elements that will be on the diagonal of the
     *              output matrix.
     * @param[out] diagonalMatrix The diagonal matrix created by this function.
     *
     * @post diagonalMatrix.distributed() == v.distributed()
     * @post diagonalMatrix.numRows() == v.dim()
     */
    Matrix IdentityMatrixFactory(const Vector &v);

}

#endif
