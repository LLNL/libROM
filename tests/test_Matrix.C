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

// Description: This source file is a test runner that uses the Google Test
// Framework to run unit tests.

#include <iostream>

#ifdef CAROM_HAS_GTEST
#include<gtest/gtest.h>
#include <mpi.h>
#include "../Matrix.h"

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
  SUCCEED();
}

/**
 *  Test methods that do not require assigning data to the matrix. These
 *  methods are:
 *
 *  TODO(oxberry1@llnl.gov): Do more exhaustive testing to test all branches
 *  of each method. For now, simple tests are better than nothing.
 *
 *  * Matrix(int, int, bool)
 *
 *  * bool distributed() const
 *
 *  * bool balanced() const
 *
 *  * int numRows() const
 *
 *  * int numCols() const
 *
 *  * void setSize(int, int)
 *
 */
TEST(MatrixSerialTest, Test_distributed)
{
  CAROM::Matrix f(2, 2, false);
  EXPECT_FALSE(f.distributed());

  CAROM::Matrix t(2, 2, true);
  EXPECT_TRUE(t.distributed());
}

TEST(MatrixSerialTest, Test_balanced)
{
  CAROM::Matrix two_by_two(2, 2, false);
  EXPECT_TRUE(two_by_two.balanced());
}

TEST(MatrixSerialTest, Test_numRows)
{
  CAROM::Matrix two_by_two(2, 2, false);
  CAROM::Matrix two_by_three(2, 3, false);
  CAROM::Matrix three_by_two(3, 2, false);

  EXPECT_EQ(two_by_two.numRows()  , 2);
  EXPECT_EQ(two_by_three.numRows(), 2);
  EXPECT_EQ(three_by_two.numRows(), 3);
}

TEST(MatrixSerialTest, Test_numColumns)
{
  CAROM::Matrix two_by_two(2, 2, false);
  CAROM::Matrix two_by_three(2, 3, false);
  CAROM::Matrix three_by_two(3, 2, false);

  EXPECT_EQ(two_by_two.numColumns()  , 2);
  EXPECT_EQ(two_by_three.numColumns(), 3);
  EXPECT_EQ(three_by_two.numColumns(), 2);
}

TEST(MatrixSerialTest, Test_setSize)
{
  CAROM::Matrix one_by_one(1, 1, false);

  EXPECT_EQ(one_by_one.numRows(), 1);
  EXPECT_EQ(one_by_one.numColumns(), 1);

  one_by_one.setSize(2, 2);
  EXPECT_EQ(one_by_one.numRows(), 2);
  EXPECT_EQ(one_by_one.numColumns(), 2);
}


/** Test methods that require assigning data to the matrix
 *
 *  TODO(oxberry1@llnl.gov): Do more exhaustive testing to test all branches
 *  of each method. For now, simple tests are better than nothing.
 *
 *  * const double& operator() (int, int) const
 *
 *  * double& operator() (int, int)
 *
 *  * const double& item(int, int) const
 *
 *  * double& item(int, int)
 *
 *  * Matrix(double*, int, int, bool, bool)
 *
 *  * Matrix(const Matrix& other)
 *
 *  * Matrix& operator= (const Matrix&)
 *
 *  * Matrix* mult(const Matrix&) const
 *
 *  * Matrix* mult(const Matrix*) const
 *
 *  * void mult(const Matrix&, Matrix*&) const
 *
 *  * void mult(const Matrix&, Matrix&) const
 *
 *  * Matrix* transposeMult(const Matrix&) const
 *
 *  * Matrix* transposeMult(const Matrix*) const
 *
 *  * void transposeMult(const Matrix&, Matrix*&) const
 *
 *  * void transposeMult(const Matrix&, Matrix&) const
 *
 *  * Matrix* inverse() const
 *
 *  * void inverse(Matrix*&) const
 *
 *  * void inverse(Matrix&) const
 *
 *  * void inverse()
 *
 *  * void qrcp_pivots_transpose(int*, int*, int) const
 */

TEST(MatrixSerialTest, Test_5arg_constructor_const_function_call)
{
  /**
   *  Build matrix [-2.0   1.0]
   *               [ 1.0  -2.0]
   *
   */
  double symmetric[4] = {-2.0, 1.0, 1.0, -2.0};
  const CAROM::Matrix symmetric_matrix(symmetric, 2, 2, false, true);

  EXPECT_DOUBLE_EQ(symmetric_matrix(0, 0), -2.0);
  EXPECT_DOUBLE_EQ(symmetric_matrix(0, 1),  1.0);
  EXPECT_DOUBLE_EQ(symmetric_matrix(1, 0),  1.0);
  EXPECT_DOUBLE_EQ(symmetric_matrix(1, 1), -2.0);

  /**
   *  Build matrix [ 1.0   0.0]
   *               [ 1.0   1.0]
   *
   */
  double asymmetric[4] = {1.0, 1.0, 0.0, 1.0};
  const CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

  EXPECT_DOUBLE_EQ(asymmetric_matrix(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix(0, 1), 1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix(1, 1), 1.0);
}

TEST(MatrixSerialTest, Test_5arg_constructor_nonconst_function_call)
{
  /**
   *  Build matrix [-2.0   1.0]
   *               [ 1.0  -2.0]
   *
   */
  double symmetric[4] = {-2.0, 1.0, 1.0, -2.0};
  CAROM::Matrix symmetric_matrix(symmetric, 2, 2, false, true);

  EXPECT_DOUBLE_EQ(symmetric_matrix(0, 0), -2.0);
  EXPECT_DOUBLE_EQ(symmetric_matrix(0, 1),  1.0);
  EXPECT_DOUBLE_EQ(symmetric_matrix(1, 0),  1.0);
  EXPECT_DOUBLE_EQ(symmetric_matrix(1, 1), -2.0);

  /**
   *  Build matrix [ 1.0   0.0]
   *               [ 1.0   1.0]
   *
   */
  double asymmetric[4] = {1.0, 1.0, 0.0, 1.0};
  CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

  EXPECT_DOUBLE_EQ(asymmetric_matrix(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix(0, 1), 1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix(1, 1), 1.0);
}

TEST(MatrixSerialTest, Test_5arg_constructor_const_item)
{
  /**
   *  Build matrix [-2.0   1.0]
   *               [ 1.0  -2.0]
   *
   */
  double symmetric[4] = {-2.0, 1.0, 1.0, -2.0};
  const CAROM::Matrix symmetric_matrix(symmetric, 2, 2, false, true);

  EXPECT_DOUBLE_EQ(symmetric_matrix.item(0, 0), -2.0);
  EXPECT_DOUBLE_EQ(symmetric_matrix.item(0, 1),  1.0);
  EXPECT_DOUBLE_EQ(symmetric_matrix.item(1, 0),  1.0);
  EXPECT_DOUBLE_EQ(symmetric_matrix.item(1, 1), -2.0);

  /**
   *  Build matrix [ 1.0   0.0]
   *               [ 1.0   1.0]
   *
   */
  double asymmetric[4] = {1.0, 1.0, 0.0, 1.0};
  const CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

  EXPECT_DOUBLE_EQ(asymmetric_matrix.item(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix.item(0, 1), 1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix.item(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix.item(1, 1), 1.0);
}

TEST(MatrixSerialTest, Test_5arg_constructor_nonconst_item)
{
  /**
   *  Build matrix [-2.0   1.0]
   *               [ 1.0  -2.0]
   *
   */
  double symmetric[4] = {-2.0, 1.0, 1.0, -2.0};
  CAROM::Matrix symmetric_matrix(symmetric, 2, 2, false, true);

  EXPECT_DOUBLE_EQ(symmetric_matrix.item(0, 0), -2.0);
  EXPECT_DOUBLE_EQ(symmetric_matrix.item(0, 1),  1.0);
  EXPECT_DOUBLE_EQ(symmetric_matrix.item(1, 0),  1.0);
  EXPECT_DOUBLE_EQ(symmetric_matrix.item(1, 1), -2.0);

  /**
   *  Build matrix [ 1.0   0.0]
   *               [ 1.0   1.0]
   *
   */
  double asymmetric[4] = {1.0, 1.0, 0.0, 1.0};
  CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

  EXPECT_DOUBLE_EQ(asymmetric_matrix.item(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix.item(0, 1), 1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix.item(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix.item(1, 1), 1.0);
}

TEST(MatrixSerialTest, Test_copy_constructor)
{
  /**
   *  Build matrix [ 1.0   0.0]
   *               [ 1.0   1.0]
   *
   */
  double asymmetric[4] = {1.0, 1.0, 0.0, 1.0};
  CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);
  CAROM::Matrix asymmetric_matrix2(asymmetric_matrix);

  EXPECT_EQ(asymmetric_matrix2.numRows(), 2);
  EXPECT_EQ(asymmetric_matrix2.numColumns(), 2);
  EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(0, 1), 1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(1, 1), 1.0);
}

TEST(MatrixSerialTest, Test_copy_assignment_operator)
{
  /**
   *  Build matrix [ 1.0   0.0]
   *               [ 1.0   1.0]
   *
   */
  double asymmetric[4] = {1.0, 1.0, 0.0, 1.0};
  CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);
  CAROM::Matrix asymmetric_matrix2 = asymmetric_matrix;

  EXPECT_EQ(asymmetric_matrix2.numRows(), 2);
  EXPECT_EQ(asymmetric_matrix2.numColumns(), 2);
  EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(0, 1), 1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(1, 1), 1.0);
}

TEST(MatrixSerialTest, Test_assignment_operator)
{
  /**
   *  Build matrix [ 1.0   0.0]
   *               [ 1.0   1.0]
   *
   */
  double asymmetric[4] = {1.0, 1.0, 0.0, 1.0};
  CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);
  CAROM::Matrix asymmetric_matrix2(2, 2, false);
  asymmetric_matrix2 = asymmetric_matrix;

  EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(0, 1), 1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(1, 1), 1.0);
}

TEST(MatrixSerialTest, Test_pMatrix_mult_reference)
{
  /**
   *  Build matrix [ 1.0   0.0]
   *               [ 1.0   1.0]
   *
   */
  double asymmetric[4] = {1.0, 1.0, 0.0, 1.0};
  const CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);
  CAROM::Matrix asymmetric_matrix2(asymmetric, 2, 2, false, true);
  CAROM::Matrix *result = NULL;

  /**
   *  [ 1.0   0.0]  *  [ 1.0   0.0]  =  [1.0   0.0]
   *  [ 1.0   1.0]     [ 1.0   1.0]     [2.0   1.0]
   *
   */
  result = asymmetric_matrix.mult(asymmetric_matrix2);
  EXPECT_EQ(result->numRows(), 2);
  EXPECT_EQ(result->numColumns(), 2);
  EXPECT_DOUBLE_EQ(result->item(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(result->item(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(result->item(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(result->item(1, 1), 1.0);

  delete result;
}

TEST(MatrixSerialTest, Test_pMatrix_mult_pointer)
{
  /**
   *  Build matrix [ 1.0   0.0]
   *               [ 1.0   1.0]
   *
   */
  double asymmetric[4] = {1.0, 1.0, 0.0, 1.0};
  const CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);
  CAROM::Matrix asymmetric_matrix2(asymmetric, 2, 2, false, true);
  CAROM::Matrix *result = NULL;

  /**
   *  [ 1.0   0.0]  *  [ 1.0   0.0]  =  [1.0   0.0]
   *  [ 1.0   1.0]     [ 1.0   1.0]     [2.0   1.0]
   *
   */
  result = asymmetric_matrix.mult(&asymmetric_matrix2);
  EXPECT_EQ(result->numRows(), 2);
  EXPECT_EQ(result->numColumns(), 2);
  EXPECT_DOUBLE_EQ(result->item(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(result->item(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(result->item(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(result->item(1, 1), 1.0);

  delete result;
}

TEST(MatrixSerialTest, Test_void_mult_output_reference)
{
  /**
   *  Build matrix [ 1.0   0.0]
   *               [ 1.0   1.0]
   *
   */
  double asymmetric[4] = {1.0, 1.0, 0.0, 1.0};
  const CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);
  CAROM::Matrix asymmetric_matrix2(asymmetric, 2, 2, false, true);
  CAROM::Matrix result(2, 2, false);

  /**
   *  [ 1.0   0.0]  *  [ 1.0   0.0]  =  [1.0   0.0]
   *  [ 1.0   1.0]     [ 1.0   1.0]     [2.0   1.0]
   *
   */
  asymmetric_matrix.mult(asymmetric_matrix, result);
  EXPECT_EQ(result.numRows(), 2);
  EXPECT_EQ(result.numColumns(), 2);
  EXPECT_DOUBLE_EQ(result.item(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(result.item(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(result.item(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(result.item(1, 1), 1.0);
}

TEST(MatrixSerialTest, Test_void_inverse_reference)
{
 /**
   *  Build matrix [ 1.0   0.0]
   *               [ 1.0   1.0]
   *
   */
  double asymmetric[4] = {1.0, 1.0, 0.0, 1.0};
  CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);
  CAROM::Matrix asymmetric_matrix_inverse(2, 2, false);

  asymmetric_matrix.inverse(asymmetric_matrix_inverse);
  EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse.item(0, 0),  1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse.item(0, 1), -1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse.item(1, 0),  0.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse.item(1, 1),  1.0);
}

TEST(MatrixSerialTest, Test_void_inverse_pointer_reference)
{
 /**
   *  Build matrix [ 1.0   0.0]
   *               [ 1.0   1.0]
   *
   */
  double asymmetric[4] = {1.0, 1.0, 0.0, 1.0};
  CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

  /**
   *  [ 1.0   0.0] ^ (-1)   =   [ 1.0   0.0]
   *  [ 1.0   1.0]              [-1.0   1.0]
   *
   */
  CAROM::Matrix *asymmetric_matrix_inverse = NULL;
  asymmetric_matrix_inverse = new CAROM::Matrix(2, 2, false);
  asymmetric_matrix.inverse(asymmetric_matrix_inverse);

  EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse->item(0, 0),  1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse->item(0, 1), -1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse->item(1, 0),  0.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse->item(1, 1),  1.0);

  delete asymmetric_matrix_inverse;
}

TEST(MatrixSerialTest, Test_void_inverse_in_place)
{
  /**
   *  Build matrix [ 1.0   0.0]
   *               [ 1.0   1.0]
   *
   */
  double asymmetric[4] = {1.0, 1.0, 0.0, 1.0};
  CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

  /**
   *  [ 1.0   0.0] ^ (-1)   =   [ 1.0   0.0]
   *  [ 1.0   1.0]              [-1.0   1.0]
   *
   */
  asymmetric_matrix.inverse();

  EXPECT_DOUBLE_EQ(asymmetric_matrix.item(0, 0),  1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix.item(0, 1), -1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix.item(1, 0),  0.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix.item(1, 1),  1.0);
}

TEST(MatrixSerialTest, Test_pMatrix_inverse)
{
  /**
   *  Build matrix [ 1.0   0.0]
   *               [ 1.0   1.0]
   *
   */
  double asymmetric[4] = {1.0, 1.0, 0.0, 1.0};
  const CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

  /**
   *  [ 1.0   0.0] ^ (-1)   =   [ 1.0   0.0]
   *  [ 1.0   1.0]              [-1.0   1.0]
   *
   */
  CAROM::Matrix* asymmetric_matrix_inverse = NULL;
  asymmetric_matrix_inverse = asymmetric_matrix.inverse();

  EXPECT_EQ(asymmetric_matrix_inverse->numRows(), 2);
  EXPECT_EQ(asymmetric_matrix_inverse->numColumns(), 2);
  EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse->item(0, 0),  1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse->item(0, 1), -1.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse->item(1, 0),  0.0);
  EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse->item(1, 1),  1.0);

  delete asymmetric_matrix_inverse;
}

/* Use second difference matrix as one fake Matrix for testing */
class SecondDifferenceMatrix : public CAROM::Matrix
{
 public:

  /**
   *  @brief Constructor.
   *
   */
  SecondDifferenceMatrix
    (int size)
    : Matrix(size, size, false)
    {
      /** A "second difference" matrix is a tridiagonal matrix that looks like:
       *
       *  [ 2 -1         ]
       *  [-1  2 -1      ]
       *  [   .. .. ..   ]
       *  [      -1  2 -1]
       *  [         -1  2]
       *
       *  which looks like the stencil from second-order central
       *  differencing (or the stiffness matrix for finite elements on a
       *  subspace of H1 with linear basis functions) for Poisson's
       *  equation, modulo boundary conditions (that would be applied to
       *  the first and last rows).
       */

      /** For degenerate case where size == 1, matrix should be one
	  entry: [2]
       */
      item(0, 0) = 2;
      if (size == 1) { return; }

      /** Otherwise, build rest of the first row. */
      item(0, 1) = -1;

      /** Build middle rows of second difference matrix. If size <= 2, this
       *  loop is not executed.
       */
      for (int row = 1; row < size - 1; row++)
      {
	item(row, row - 1) = -1;
	item(row, row    ) =  2;
	item(row, row + 1) = -1;
      }

      /** Build last row of second difference matrix. If size == 1, these
       *  statements are not executed.
       */
      item(size - 1, size - 2) = -1;
      item(size - 1, size - 1) =  2;
      return;
    }

  /**
   * @brief Destructor.
   */
  ~SecondDifferenceMatrix
    ()
    {
    }
};

// Test fixture for testing
TEST(MatrixSerialTest, Test_qrcp_pivots_transpose)
{
  // Allocate space for second_difference matrix
  const int size = 4;
  SecondDifferenceMatrix second_difference(size);

  // Allocate memory for QRCP arguments
  int* row_pivot                 =       new int[size];
  int* row_pivot_owner           =       new int[size];
  const int row_pivots_requested =               size ;

  // Compute pivots
  second_difference.qrcp_pivots_transpose(row_pivot,
					  row_pivot_owner,
					  row_pivots_requested);

  // All row pivot owners should be on the current rank;
  // all row pivots shouold be less than the size of the matrix
  int is_mpi_initialized, is_mpi_finalized;
  CAROM_ASSERT(MPI_Initialized(&is_mpi_initialized) == MPI_SUCCESS);
  CAROM_ASSERT(MPI_Finalized(&is_mpi_finalized) == MPI_SUCCESS);
  int my_rank = 0;
  if(is_mpi_initialized && !is_mpi_finalized) {
    const MPI_Comm my_comm = MPI_COMM_WORLD;
    CAROM_ASSERT(MPI_Comm_rank(my_comm, &my_rank) == MPI_SUCCESS);
  }

  for (int i = 0; i < row_pivots_requested; i++) {
    EXPECT_EQ(row_pivot_owner[i], my_rank);
    EXPECT_TRUE(row_pivot[i] < size);
  }

  // Row pivots should be known values
  EXPECT_EQ(row_pivot[0], 1);
  EXPECT_EQ(row_pivot[1], 3);
  EXPECT_EQ(row_pivot[2], 0);
  EXPECT_EQ(row_pivot[3], 2);

  // Free allocated arrays
  delete [] row_pivot;
  delete [] row_pivot_owner;
}

/**
 * Test methods that take CAROM::Vector objects as input and/or output
 *
 * * Vector* mult(const Vector&)
 *
 * * Vector* mult(const Vector*)
 *
 * * void mult(const Vector&, Vector*&)
 *
 * * void mult(const Vector&, Vector&) const
 *
 * * void multPlus(Vector&, const Vector&, double) const
 *
 * * Vector* transposeMult(const Vector&) const
 *
 * * Vector* transposeMult(const Vector*) const
 *
 * * void transposeMult(const Vector&, Vector*&) const
 *
 * * void transposeMult(const Vector&, Vector&) const
 *
 */

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#else // #ifndef CAROM_HAS_GTEST
int main()
{
  std::cout << "libROM was compiled without Google Test support, so unit "
	    << "tests have been disabled. To enable unit tests, compile "
	    << "libROM with Google Test support." << std::endl;
}
#endif // #endif CAROM_HAS_GTEST
