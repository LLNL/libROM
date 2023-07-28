/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
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
#include <vector>
#include <complex>
#include <string>

namespace CAROM {

/**
 * Class Matrix is a simple matrix class in which the rows may be distributed
 * across multiple processes. This class supports only the basic operations that
 * are needed by the SVD library.
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
     * @param[in] randomized If true the matrix will be a standard normally
                             distributed random matrix.
     */
    Matrix(
        int num_rows,
        int num_cols,
        bool distributed,
        bool randomized = false);

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
     * @param[in] copy_data If true the matrix allocates its own storage and
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
     * reallocates storage if needed. All values are initialized to zero.
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

            // Allocate new array and initialize all values to zero.
            d_mat = new double [new_size] {0.0};
            d_alloc_size = new_size;
        }
        d_num_rows = num_rows;
        d_num_cols = num_cols;
        if (d_distributed) {
            calculateNumDistributedRows();
        }
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
     * @brief Returns the number of rows of the Matrix across all processors.
     *
     * @return The number of rows of the Matrix across all processors.
     */
    int
    numDistributedRows() const
    {
        if (!d_distributed) {
            return d_num_rows;
        }
        return d_num_distributed_rows;
    }

    /**
     * @brief Returns the number of columns in the Matrix. This method will
      *       return the same value from each processor.
     *
     * @return The number of columns of the Matrix.
     */
    int
    numColumns() const
    {
        return d_num_cols;
    }

    /**
     * @brief Get the first N columns of a matrix.
     *
     * @pre 0 < n < numColumns()
     *
     * @param[in] n The number of columns to return.
     *
     * @return The truncated Matrix.
     */
    Matrix*
    getFirstNColumns(int n) const;

    /**
     * @brief Get the first N columns of a matrix.
     *
     * @pre result.distributed() == distributed()
     * @pre 0 < n < numColumns()
     *
     * @param[in] n The number of columns to return.
     *
     * @param[out] result The truncated Matrix.
     */
    void
    getFirstNColumns(
        int n,
        Matrix*& result) const;

    /**
     * @brief Get the first N columns of a matrix.
     *
     * @pre result.distributed() == distributed()
     * @pre 0 < n < numColumns()
     *
     * @param[in] n The number of columns to return.
     *
     * @param[out] result The truncated Matrix.
     */
    void
    getFirstNColumns(
        int n,
        Matrix& result) const;

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
        CAROM_VERIFY(other != 0);
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
        CAROM_VERIFY(other != 0);
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
    * @param[in] this_row The row of the matrix to multiple with other.
    * @param[in] other The Vector to multiply with this_row of the matrix.
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
    * @param[in] this_row The row of the matrix to multiple with other.
    * @param[in] other The Vector to multiply with this_row of the matrix.
    * @param[out] other The product Vector.
    */
    void
    pointwise_mult(
        int this_row,
        Vector& other) const;

    /**
     * @brief Multiplies two matrices element-wise.
     *
     * @pre result.distributed() == distributed()
     * @pre distributed() == other.distributed()
     * @pre numRows() == other.numRows()
     * @pre numColumns() == other.numColumns()
     *
     * @param[in] other The Matrix to multiply with this.
     *
     * @return The product Matrix.
     */
    Matrix*
    elementwise_mult(
        const Matrix& other) const
    {
        Matrix* result = 0;
        elementwise_mult(other, result);
        return result;
    }

    /**
     * @brief Multiplies two matrices element-wise and returns the product,
     * pointer version.
     *
     * @pre other != 0
     * @pre distributed() == other.distributed()
     * @pre numRows() == other.numRows()
     * @pre numColumns() == other.numColumns()
     *
     * @param[in] other The Matrix to multiply with this.
     *
     * @return The product Matrix.
     */
    Matrix*
    elementwise_mult(
        const Matrix* other) const
    {
        CAROM_VERIFY(other != 0);
        return elementwise_mult(*other);
    }

    /**
     * @brief Multiplies two matrices element-wise and fills result with the
     * answer.
     *
     * @pre result == 0 || result->distributed() == distributed()
     * @pre distributed() == other.distributed()
     * @pre numRows() == other.numRows()
     * @pre numColumns() == other.numColumns()
     *
     * @param[in] other The Matrix to multiply with this.
     * @param[out] result The product Matrix.
     */
    void
    elementwise_mult(
        const Matrix& other,
        Matrix*& result) const;

    /**
     * @brief Multiplies two matrices element-wise and fills result with the
     * answer.
     *
     * @pre result == 0 || result->distributed() == distributed()
     * @pre distributed() == other.distributed()
     * @pre numRows() == other.numRows()
     * @pre numColumns() == other.numColumns()
     *
     * @param[in] other The Matrix to multiply with this.
     * @param[out] result The product Matrix.
     */
    void
    elementwise_mult(
        const Matrix& other,
        Matrix& result) const;

    /**
     * @brief Square every element in the matrix.
     *
     * @return The product Matrix.
     */
    Matrix*
    elementwise_square() const
    {
        Matrix* result = 0;
        elementwise_square(result);
        return result;
    }

    /**
     * @brief Square every element in the matrix.
     *
     * @pre result == 0 || result->distributed() == distributed()
     * @param[out] result The product Matrix.
     */
    void
    elementwise_square(
        Matrix*& result) const;

    /**
     * @brief Square every element in the matrix.
     *
     * @pre result == 0 || result->distributed() == distributed()
     * @param[out] result The product Matrix.
     */
    void
    elementwise_square(
        Matrix& result) const;

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
        CAROM_VERIFY(other != 0);
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
        CAROM_VERIFY(other != 0);
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
     * @brief Returns a column of the matrix (not owned by Matrix).
     *
     * @return A column of the matrix (not owned by Matrix).
     */
    Vector*
    getColumn(int column) const
    {
        Vector* result = 0;
        getColumn(column, result);
        return result;
    }

    /**
     * @brief Returns a column of the matrix (not owned by Matrix).
     *
     * @return A column of the matrix (not owned by Matrix).
     */
    void
    getColumn(int column,
              Vector*& result) const;

    /**
     * @brief Returns a column of the matrix (not owned by Matrix).
     *
     * @return A column of the matrix (not owned by Matrix).
     */
    void
    getColumn(int column,
              Vector& result) const;

    /**
     * @brief Replaces this Matrix with its transpose (in place),
     * in the serial square case only.
     *
     * @pre !distributed()
     * @pre numRows() == numColumns()
     */
    void transpose();

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

    /**
     * @brief Computes and returns the Q of the QR factorization of this.
     *
     * @return The Q of the QR factorization of this.
     */
    Matrix*
    qr_factorize() const;

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
     * @param[in] pivots_requested The number of pivots requested, must be less
     *                             than or equal to the number of rows of this
     *                             Matrix.
     */
    void
    qrcp_pivots_transpose(int* row_pivot,
                          int* row_pivot_owner,
                          int  pivots_requested) const;

    /**
     * @brief Orthogonalizes the matrix.
     */
    void
    orthogonalize();

    /**
     * @brief Rescale every matrix row by its maximum absolute value.
     */
    void
    rescale_rows_max();

    /**
     * @brief Rescale every matrix column by its maximum absolute value.
     */
    void
    rescale_cols_max();

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
    void print(const char * prefix) const;

    /**
     * @brief write Matrix into (a) HDF file(s).
     *
     * @param[in] base_file_name The base part of the file name.
     *
     */
    void write(const std::string& base_file_name) const;

    /**
     * @brief read Matrix into (a) HDF file(s).
     *
     * @param[in] base_file_name The base part of the file name.
     *
     */
    void read(const std::string& base_file_name);

    /**
     * @brief read a single rank of a distributed Matrix into (a) HDF file(s).
     *
     * @param[in] base_file_name The base part of the file name.
     * @param[in] rank The rank to read from.
     *
     */
    void local_read(const std::string& base_file_name, int rank);

    /**
     * @brief Get the matrix data as a pointer.
     */
    double *getData() const
    {
        return d_mat;
    }

    /**
     * @brief Distribute this matrix rows among MPI processes,
     * based on the specified local number of rows.
     * This becomes distributed after this function is executed.
     *
     * @pre !distributed()
     * @pre d_owns_data
     *
     * @param[in] local_num_rows number of rows for local MPI rank.
     */
    void distribute(const int &local_num_rows);

    /**
     * @brief Gather all the distributed rows among MPI processes.
     * This becomes not distributed after this function is executed.
     *
     * @pre distributed()
     * @pre d_owns_data
     */
    void gather();

private:
    /**
     * @brief Compute number of rows across all processors.
     *
     * @pre distributed()
     */
    void
    calculateNumDistributedRows();

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
     * @brief The storage for the Matrix's values on this processor.
     */
    double* d_mat;

    /**
     * @brief The rows in the Matrix that are on this processor.
     */
    int d_num_rows;

    /**
     * @brief The rows in the Matrix across all processors.
     */
    int d_num_distributed_rows;

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
     * @brief The current MPI rank. If MPI is not initialized, equal to 0.
     */
    int d_rank;

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
 *
 * @post diagonalMatrix.distributed() == v.distributed()
 * @post diagonalMatrix.numRows() == v.dim()
 *
 * @return The diagonal matrix created by this function.
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
 *
 * @post diagonalMatrix.distributed() == v.distributed()
 * @post diagonalMatrix.numRows() == v.dim()
 *
 * @return The identity matrix created by this function.
 */
Matrix IdentityMatrixFactory(const Vector &v);

/**
 * struct EigenPair is a struct to hold the real eigenvectors
 * of a matrix along with its real eigenvalues.
 */
struct EigenPair {

    /**
     * @brief The eigenvectors in matrix form.
     */
    Matrix* ev;

    /**
     * @brief The corresponding eigenvalues.
     */
    std::vector<double> eigs;
};

/**
 * struct ComplexEigenPair is a struct to hold the real and imaginary eigenvectors
 * of a matrix along with its eigenvalues.
 */
struct ComplexEigenPair {

    /**
     * @brief The real component of the eigenvectors in matrix form.
     */
    Matrix* ev_real;

    /**
     * @brief The imaginary component of the eigenvectors in matrix form.
     */
    Matrix* ev_imaginary;

    /**
     * @brief The corresponding complex eigenvalues.
     */
    std::vector<std::complex<double>> eigs;
};

/**
 * struct ComplexEigenPair is a struct to hold the singular value decomposition
 * of an undistributed matrix.
 */
struct SerialSVDDecomposition
{
    /**
     * @brief The left singular vectors.
     */
    Matrix* U;

    /**
     * @brief The singular values.
     */
    Vector* S;

    /**
     * @brief The right singular vectors.
     */
    Matrix* V;
};

/**
 * @brief Computes the SVD of a undistributed matrix in serial.
 *
 * @param[in] A The MxN undistributed matrix to be eigendecomposed.
 * @param[in] U The matrix to write the left singular vectors into.
 * @param[in] S The vector to write the singular values into.
 * @param[in] V The matrix to write the right singular vectors into.
 *
 * @return The singular value decomposition within the input parameters U, S,
 *         and V.
 */
void SerialSVD(Matrix* A,
               Matrix* U,
               Vector* S,
               Matrix* V);

/**
 * @brief Computes the SVD of a undistributed matrix in serial.
 *
 * @param[in] A The MxN undistributed matrix to be eigendecomposed.
 *
 * @return The singular value decomposition within a struct. The matrices
 *         and vector contained within the returning struct must be destroyed by
 *         the user.
 */
struct SerialSVDDecomposition SerialSVD(Matrix* A);

/**
 * @brief Computes the eigenvectors/eigenvalues of an NxN real symmetric matrix.
 *
 * @param[in] A The NxN real symmetric matrix to be eigendecomposed.
 *
 * @return The eigenvectors and eigenvalues of the eigensolve. The eigenvector matrices
 *         contained within the returning struct must be destroyed by the user.
 */
struct EigenPair SymmetricRightEigenSolve(Matrix* A);

/**
 * @brief Computes the eigenvectors/eigenvalues of an NxN real nonsymmetric matrix.
 *
 * @param[in] A The NxN real nonsymmetric matrix to be eigendecomposed.
 *
 * @return The eigenvectors and eigenvalues of the eigensolve. The eigenvector matrices
 *         contained within the returning struct must be destroyed by the user.
 */
struct ComplexEigenPair NonSymmetricRightEigenSolve(Matrix* A);

Matrix* SpaceTimeProduct(const CAROM::Matrix* As, const CAROM::Matrix* At,
                         const CAROM::Matrix* Bs, const CAROM::Matrix* Bt,
                         const std::vector<double> *tscale=NULL,
                         const bool At0at0=false, const bool Bt0at0=false,
                         const bool lagB=false, const bool skip0=false);
}

#endif
