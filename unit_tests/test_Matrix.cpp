/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: This source file is a test runner that uses the Google Test
// Framework to run unit tests on the CAROM::Matrix class.

#include <iostream>

#ifdef CAROM_HAS_GTEST
#include<gtest/gtest.h>
#include <mpi.h>
#include "linalg/Matrix.h"
#include "utils/mpi_utils.h"

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

    EXPECT_EQ(two_by_two.numRows(), 2);
    EXPECT_EQ(two_by_three.numRows(), 2);
    EXPECT_EQ(three_by_two.numRows(), 3);
}

TEST(MatrixSerialTest, Test_numColumns)
{
    CAROM::Matrix two_by_two(2, 2, false);
    CAROM::Matrix two_by_three(2, 3, false);
    CAROM::Matrix three_by_two(3, 2, false);

    EXPECT_EQ(two_by_two.numColumns(), 2);
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
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

    EXPECT_DOUBLE_EQ(asymmetric_matrix.item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix.item(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix.item(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix.item(1, 1), 1.0);
}

TEST(MatrixSerialTest, Test_copy_constructor)
{
    /**
     *  Build matrix [ 1.0   0.0]
     *               [ 1.0   1.0]
     *
     */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);
    CAROM::Matrix asymmetric_matrix2(asymmetric_matrix);

    EXPECT_EQ(asymmetric_matrix2.numRows(), 2);
    EXPECT_EQ(asymmetric_matrix2.numColumns(), 2);
    EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(1, 1), 1.0);
}

TEST(MatrixSerialTest, Test_copy_assignment_operator)
{
    /**
     *  Build matrix [ 1.0   0.0]
     *               [ 1.0   1.0]
     *
     */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);
    CAROM::Matrix asymmetric_matrix2 = asymmetric_matrix;

    EXPECT_EQ(asymmetric_matrix2.numRows(), 2);
    EXPECT_EQ(asymmetric_matrix2.numColumns(), 2);
    EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(1, 1), 1.0);
}

TEST(MatrixSerialTest, Test_assignment_operator)
{
    /**
     *  Build matrix [ 1.0   0.0]
     *               [ 1.0   1.0]
     *
     */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);
    CAROM::Matrix asymmetric_matrix2(2, 2, false);
    asymmetric_matrix2 = asymmetric_matrix;

    EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix2.item(1, 1), 1.0);
}

TEST(MatrixSerialTest, Test_get_first_n_columns)
{
    /**
     *  Build matrix [ 0.0   1.0   2.0   3.0]
     *               [ 4.0   5.0   6.0   7.0]
     *               [ 8.0   9.0  10.0  11.0]
     *               [12.0  13.0  14.0  15.0]
     */
    double d_mat[16] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                        11.0, 12.0, 13.0, 14.0, 15.0
                       };
    CAROM::Matrix matrix(d_mat, 4, 4, false);
    CAROM::Matrix* truncated_matrix = matrix.getFirstNColumns(2);

    EXPECT_EQ(truncated_matrix->numRows(), 4);
    EXPECT_EQ(truncated_matrix->numColumns(), 2);

    EXPECT_DOUBLE_EQ(truncated_matrix->item(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(truncated_matrix->item(0, 1), 1.0);
    EXPECT_DOUBLE_EQ(truncated_matrix->item(1, 0), 4.0);
    EXPECT_DOUBLE_EQ(truncated_matrix->item(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(truncated_matrix->item(2, 0), 8.0);
    EXPECT_DOUBLE_EQ(truncated_matrix->item(2, 1), 9.0);
    EXPECT_DOUBLE_EQ(truncated_matrix->item(3, 0), 12.0);
    EXPECT_DOUBLE_EQ(truncated_matrix->item(3, 1), 13.0);
}

TEST(MatrixSerialTest, Test_pMatrix_mult_reference)
{
    /**
     *  Build matrix [ 1.0   0.0]
     *               [ 1.0   1.0]
     *
     */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
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
    EXPECT_DOUBLE_EQ(result->item(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(result->item(1, 0), 2.0);
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
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
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
    EXPECT_DOUBLE_EQ(result->item(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(result->item(1, 0), 2.0);
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
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    const CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);
    CAROM::Matrix asymmetric_matrix2(asymmetric, 2, 2, false, true);
    CAROM::Matrix result(2, 2, false);

    /**
     *  [ 1.0   0.0]  *  [ 1.0   0.0]  =  [1.0   0.0]
     *  [ 1.0   1.0]     [ 1.0   1.0]     [2.0   1.0]
     *
     */
    asymmetric_matrix.mult(asymmetric_matrix2, result);
    EXPECT_EQ(result.numRows(), 2);
    EXPECT_EQ(result.numColumns(), 2);
    EXPECT_DOUBLE_EQ(result.item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(result.item(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(result.item(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(result.item(1, 1), 1.0);
}

TEST(MatrixSerialTest, Test_void_mult_output_pointer)
{
    /**
     *  Build matrix [ 1.0   0.0]
     *               [ 1.0   1.0]
     *
     */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    const CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);
    CAROM::Matrix asymmetric_matrix2(asymmetric, 2, 2, false, true);
    CAROM::Matrix *result = new CAROM::Matrix(2, 2, false);

    /**
     *  [ 1.0   0.0]  *  [ 1.0   0.0]  =  [1.0   0.0]
     *  [ 1.0   1.0]     [ 1.0   1.0]     [2.0   1.0]
     *
     */
    asymmetric_matrix.mult(asymmetric_matrix2, result);
    EXPECT_EQ(result->numRows(), 2);
    EXPECT_EQ(result->numColumns(), 2);
    EXPECT_DOUBLE_EQ(result->item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(result->item(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(result->item(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(result->item(1, 1), 1.0);

    delete result;
}

TEST(MatrixSerialTest, Test_pMatrix_transpose_mult_output_reference)
{
    /**
     *  Build matrix [ 1.0   0.0]
     *               [ 1.0   1.0]
     *
     */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    const CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

    /**
     *  Build identity matrix  [ 1.0   0.0]
     *                         [ 0.0   1.0]
     *
     */
    double identity[4] = {1.0, 0.0, 0.0, 1.0};
    const CAROM::Matrix identity_matrix(identity, 2, 2, false, true);
    EXPECT_DOUBLE_EQ(identity_matrix.item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(identity_matrix.item(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(identity_matrix.item(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(identity_matrix.item(1, 1), 1.0);

    /**
     *  [ 1.0   0.0]^T  *  [ 1.0   0.0]    =  [1.0   1.0]
     *  [ 1.0   1.0]       [ 0.0   1.0]       [0.0   1.0]
     *
     */
    CAROM::Matrix *result = NULL;
    result = asymmetric_matrix.transposeMult(identity_matrix);
    EXPECT_EQ(result->numRows(), 2);
    EXPECT_EQ(result->numColumns(), 2);
    EXPECT_DOUBLE_EQ(result->item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(result->item(0, 1), 1.0);
    EXPECT_DOUBLE_EQ(result->item(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(result->item(1, 1), 1.0);
    delete result;
    result = NULL;

    /**
     *  [ 1.0   0.0]^T  *  [ 1.0   0.0]    =  [1.0   0.0]
     *  [ 0.0   1.0]       [ 1.0   1.0]       [1.0   1.0]
     *
     */
    result = identity_matrix.transposeMult(asymmetric_matrix);
    EXPECT_EQ(result->numRows(), 2);
    EXPECT_EQ(result->numColumns(), 2);
    EXPECT_DOUBLE_EQ(result->item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(result->item(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(result->item(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(result->item(1, 1), 1.0);
    delete result;
    result = NULL;

    /**
     *  [ 1.0   0.0]^T  *  [ 1.0   0.0]    =  [2.0   1.0]
     *  [ 1.0   1.0]       [ 1.0   1.0]       [1.0   1.0]
     *
     */
    result = asymmetric_matrix.transposeMult(asymmetric_matrix);
    EXPECT_EQ(result->numRows(), 2);
    EXPECT_EQ(result->numColumns(), 2);
    EXPECT_DOUBLE_EQ(result->item(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(result->item(0, 1), 1.0);
    EXPECT_DOUBLE_EQ(result->item(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(result->item(1, 1), 1.0);
    delete result;
}

TEST(MatrixSerialTest, Test_pMatrix_transpose_mult_output_pointer)
{
    /**
     *  Build matrix [ 1.0   0.0]
     *               [ 1.0   1.0]
     *
     */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    const CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

    /**
     *  Build identity matrix  [ 1.0   0.0]
     *                         [ 0.0   1.0]
     *
     */
    double identity[4] = {1.0, 0.0, 0.0, 1.0};
    const CAROM::Matrix identity_matrix(identity, 2, 2, false, true);
    EXPECT_DOUBLE_EQ(identity_matrix.item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(identity_matrix.item(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(identity_matrix.item(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(identity_matrix.item(1, 1), 1.0);

    /**
     *  [ 1.0   0.0]^T  *  [ 1.0   0.0]    =  [1.0   1.0]
     *  [ 1.0   1.0]       [ 0.0   1.0]       [0.0   1.0]
     *
     */
    CAROM::Matrix *result = NULL;
    result = asymmetric_matrix.transposeMult(&identity_matrix);
    EXPECT_EQ(result->numRows(), 2);
    EXPECT_EQ(result->numColumns(), 2);
    EXPECT_DOUBLE_EQ(result->item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(result->item(0, 1), 1.0);
    EXPECT_DOUBLE_EQ(result->item(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(result->item(1, 1), 1.0);
    delete result;
    result = NULL;

    /**
     *  [ 1.0   0.0]^T  *  [ 1.0   0.0]    =  [1.0   0.0]
     *  [ 0.0   1.0]       [ 1.0   1.0]       [1.0   1.0]
     *
     */
    result = identity_matrix.transposeMult(&asymmetric_matrix);
    EXPECT_EQ(result->numRows(), 2);
    EXPECT_EQ(result->numColumns(), 2);
    EXPECT_DOUBLE_EQ(result->item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(result->item(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(result->item(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(result->item(1, 1), 1.0);
    delete result;
    result = NULL;

    /**
     *  [ 1.0   0.0]^T  *  [ 1.0   0.0]    =  [2.0   1.0]
     *  [ 1.0   1.0]       [ 1.0   1.0]       [1.0   1.0]
     *
     */
    result = asymmetric_matrix.transposeMult(&asymmetric_matrix);
    EXPECT_EQ(result->numRows(), 2);
    EXPECT_EQ(result->numColumns(), 2);
    EXPECT_DOUBLE_EQ(result->item(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(result->item(0, 1), 1.0);
    EXPECT_DOUBLE_EQ(result->item(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(result->item(1, 1), 1.0);
    delete result;
}

TEST(MatrixSerialTest, Test_void_transpose_mult_reference)
{
    /**
     *  Build matrix [ 1.0   0.0]
     *               [ 1.0   1.0]
     *
     */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    const CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

    /**
     *  Build identity matrix  [ 1.0   0.0]
     *                         [ 0.0   1.0]
     *
     */
    double identity[4] = {1.0, 0.0, 0.0, 1.0};
    const CAROM::Matrix identity_matrix(identity, 2, 2, false, true);
    EXPECT_DOUBLE_EQ(identity_matrix.item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(identity_matrix.item(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(identity_matrix.item(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(identity_matrix.item(1, 1), 1.0);

    /**
     *  [ 1.0   0.0]^T  *  [ 1.0   0.0]    =  [1.0   1.0]
     *  [ 1.0   1.0]       [ 0.0   1.0]       [0.0   1.0]
     *
     */
    CAROM::Matrix result(2, 2, false);
    asymmetric_matrix.transposeMult(identity_matrix, result);
    EXPECT_DOUBLE_EQ(result.item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(result.item(0, 1), 1.0);
    EXPECT_DOUBLE_EQ(result.item(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(result.item(1, 1), 1.0);

    /**
     *  [ 1.0   0.0]^T  *  [ 1.0   0.0]    =  [1.0   0.0]
     *  [ 0.0   1.0]       [ 1.0   1.0]       [1.0   1.0]
     *
     */
    identity_matrix.transposeMult(asymmetric_matrix, result);
    EXPECT_DOUBLE_EQ(result.item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(result.item(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(result.item(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(result.item(1, 1), 1.0);

    /**
     *  [ 1.0   0.0]^T  *  [ 1.0   0.0]    =  [2.0   1.0]
     *  [ 1.0   1.0]       [ 1.0   1.0]       [1.0   1.0]
     *
     */
    asymmetric_matrix.transposeMult(asymmetric_matrix, result);
    EXPECT_DOUBLE_EQ(result.item(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(result.item(0, 1), 1.0);
    EXPECT_DOUBLE_EQ(result.item(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(result.item(1, 1), 1.0);
}

TEST(MatrixSerialTest, Test_void_transpose_mult_pointer)
{
    /**
     *  Build matrix [ 1.0   0.0]
     *               [ 1.0   1.0]
     *
     */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    const CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

    /**
     *  Build identity matrix  [ 1.0   0.0]
     *                         [ 0.0   1.0]
     *
     */
    double identity[4] = {1.0, 0.0, 0.0, 1.0};
    const CAROM::Matrix identity_matrix(identity, 2, 2, false, true);
    EXPECT_DOUBLE_EQ(identity_matrix.item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(identity_matrix.item(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(identity_matrix.item(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(identity_matrix.item(1, 1), 1.0);

    /**
     *  [ 1.0   0.0]^T  *  [ 1.0   0.0]    =  [1.0   1.0]
     *  [ 1.0   1.0]       [ 0.0   1.0]       [0.0   1.0]
     *
     */
    CAROM::Matrix *result = NULL;
    asymmetric_matrix.transposeMult(identity_matrix, result);
    EXPECT_EQ(result->numRows(), 2);
    EXPECT_EQ(result->numColumns(), 2);
    EXPECT_DOUBLE_EQ(result->item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(result->item(0, 1), 1.0);
    EXPECT_DOUBLE_EQ(result->item(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(result->item(1, 1), 1.0);
    delete result;
    result = NULL;

    /**
     *  [ 1.0   0.0]^T  *  [ 1.0   0.0]    =  [1.0   0.0]
     *  [ 0.0   1.0]       [ 1.0   1.0]       [1.0   1.0]
     *
     */
    identity_matrix.transposeMult(asymmetric_matrix, result);
    EXPECT_EQ(result->numRows(), 2);
    EXPECT_EQ(result->numColumns(), 2);
    EXPECT_DOUBLE_EQ(result->item(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(result->item(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(result->item(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(result->item(1, 1), 1.0);
    delete result;
    result = NULL;

    /**
     *  [ 1.0   0.0]^T  *  [ 1.0   0.0]    =  [2.0   1.0]
     *  [ 1.0   1.0]       [ 1.0   1.0]       [1.0   1.0]
     *
     */
    asymmetric_matrix.transposeMult(asymmetric_matrix, result);
    EXPECT_EQ(result->numRows(), 2);
    EXPECT_EQ(result->numColumns(), 2);
    EXPECT_DOUBLE_EQ(result->item(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(result->item(0, 1), 1.0);
    EXPECT_DOUBLE_EQ(result->item(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(result->item(1, 1), 1.0);
    delete result;
}

TEST(MatrixSerialTest, Test_void_inverse_reference)
{
    /**
      *  Build matrix [ 1.0   0.0]
      *               [ 1.0   1.0]
      *
      */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);
    CAROM::Matrix asymmetric_matrix_inverse(2, 2, false);

    asymmetric_matrix.inverse(asymmetric_matrix_inverse);
    EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse.item(0, 0),  1.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse.item(0, 1),  0.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse.item(1, 0), -1.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse.item(1, 1),  1.0);
}

TEST(MatrixSerialTest, Test_void_inverse_pointer_reference)
{
    /**
      *  Build matrix [ 1.0   0.0]
      *               [ 1.0   1.0]
      *
      */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
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
    EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse->item(0, 1),  0.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse->item(1, 0), -1.0);
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
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

    /**
     *  [ 1.0   0.0] ^ (-1)   =   [ 1.0   0.0]
     *  [ 1.0   1.0]              [-1.0   1.0]
     *
     */
    asymmetric_matrix.inverse();

    EXPECT_DOUBLE_EQ(asymmetric_matrix.item(0, 0),  1.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix.item(0, 1),  0.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix.item(1, 0), -1.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix.item(1, 1),  1.0);
}

TEST(MatrixSerialTest, Test_pMatrix_inverse)
{
    /**
     *  Build matrix [ 1.0   0.0]
     *               [ 1.0   1.0]
     *
     */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
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
    EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse->item(0, 1),  0.0);
    EXPECT_DOUBLE_EQ(asymmetric_matrix_inverse->item(1, 0), -1.0);
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

        /* Zero out matrix */
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                item(i, j) = 0;
            }
        }

        /** For degenerate case where size == 1, matrix should be one
        entry: [2]
         */
        item(0, 0) = 2;
        if (size == 1) {
            return;
        }

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

/* Use permuted identity matrix as one fake Matrix for testing */
class PermutedIdentityMatrix : public CAROM::Matrix
{
public:
    /**
     * @brief Constructor
     */
    PermutedIdentityMatrix
    (int size, int* permutation, bool is_inverse=false)
        : Matrix(size, size, false)
    {
        /* Zero out matrix */
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                item(i, j) = 0;
            }
        }

        if (!is_inverse)
        {
            /* Assign ones in the columns  */
            for (int i = 0; i < size; i++)
            {
                item(i, permutation[i]) = 1;
            }
        }
        else
        {
            /* Assign ones in the rows; an inverse permutation matrix  */
            for (int i = 0; i < size; i++)
            {
                item(permutation[i], i) = 1;
            }
        }
    }

    /* TODO(oxberry1@llnl.gov): Replace with delegating c'tor once C++11 adopted*/
    PermutedIdentityMatrix
    (std::vector<int> permutation, bool is_inverse=false)
        : Matrix(static_cast<int>(permutation.size()),
                 static_cast<int>(permutation.size()),
                 false)
    {
        int size = permutation.size();

        /* Zero out matrix */
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                item(i, j) = 0;
            }
        }

        if (!is_inverse)
        {
            /* Assign ones in the columns  */
            for (int i = 0; i < size; i++)
            {
                item(i, permutation.at(i)) = 1;
            }
        }
        else
        {
            /* Assign ones in the rows; an inverse permutation matrix  */
            for (int i = 0; i < size; i++)
            {
                item(permutation.at(i), i) = 1;
            }
        }
    }

    ~PermutedIdentityMatrix
    ()
    {
    }
};

// Test fixture for testing; this test mainly tests to make sure
// that C-style indices (starting at zero) are returned, and in
// the correct order; the identity matrix should be returned
// in the same order.
TEST(MatrixSerialTest, Test_qrcp_pivots_transpose)
{
    // Allocate space for permuted_identity matrix
    const int size = 4;
    int permutation[4] = {0, 1, 2, 3};
    PermutedIdentityMatrix permuted_identity(size, permutation, false);

    // Allocate memory for QRCP arguments
    int* row_pivot                 =       new int[size];
    int* row_pivot_owner           =       new int[size];
    const int row_pivots_requested =               size ;

    // Compute pivots
    permuted_identity.qrcp_pivots_transpose(row_pivot,
                                            row_pivot_owner,
                                            row_pivots_requested);

    // All row pivot owners should be on the current rank;
    // all row pivots should be less than the size of the matrix
    int is_mpi_initialized, is_mpi_finalized;
    MPI_Initialized(&is_mpi_initialized);
    MPI_Finalized(&is_mpi_finalized);
    const MPI_Comm my_comm = MPI_COMM_WORLD;
    int my_rank = 0;
    if(is_mpi_initialized && !is_mpi_finalized) {
        MPI_Comm_rank(my_comm, &my_rank);
    }

    for (int i = 0; i < row_pivots_requested; i++) {
        EXPECT_EQ(row_pivot_owner[i], my_rank);
        EXPECT_TRUE(row_pivot[i] < size);
    }

    /* Test known row pivots */
    EXPECT_EQ(row_pivot[0], permutation[0]);
    EXPECT_EQ(row_pivot[1], permutation[1]);
    EXPECT_EQ(row_pivot[2], permutation[2]);
    EXPECT_EQ(row_pivot[3], permutation[3]);

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
TEST(MatrixSerialTest, Test_mult_Vector_reference)
{
    /**
      *  Build matrix [ 1.0   0.0]
      *               [ 1.0   1.0]
      *
      */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

    /**
     *  Build vector [ 2.0]
     *               [ 4.0]
     */
    double v_data[2] = {2.0, 4.0};
    CAROM::Vector v(v_data, 2, false, true);

    /**
      *  [ 1.0   0.0] [ 2.0] = [ 2.0]
      *  [ 1.0   1.0] [ 4.0]   [ 6.0]
      *
      */
    CAROM::Vector *w;
    w = asymmetric_matrix.mult(v);
    EXPECT_FALSE(w->distributed());
    EXPECT_EQ(w->dim(), 2);
    EXPECT_DOUBLE_EQ((*w)(0), 2.0);
    EXPECT_DOUBLE_EQ((*w)(1), 6.0);
    delete w;
}

TEST(MatrixSerialTest, Test_mult_Vector_pointer)
{
    /**
      *  Build matrix [ 1.0   0.0]
      *               [ 1.0   1.0]
      *
      */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

    /**
     *  Build vector [ 2.0]
     *               [ 4.0]
     */
    double v_data[2] = {2.0, 4.0};
    CAROM::Vector *v = new CAROM::Vector(v_data, 2, false, true);

    /**
      *  [ 1.0   0.0] [ 2.0] = [ 2.0]
      *  [ 1.0   1.0] [ 4.0]   [ 6.0]
      *
      */
    CAROM::Vector *w;
    w = asymmetric_matrix.mult(v);
    EXPECT_FALSE(w->distributed());
    EXPECT_EQ(w->dim(), 2);
    EXPECT_DOUBLE_EQ((*w)(0), 2.0);
    EXPECT_DOUBLE_EQ((*w)(1), 6.0);
    delete w;
    delete v;
}

TEST(MatrixSerialTest, Test_mult_Vector_Vector_pointer)
{
    /**
      *  Build matrix [ 1.0   0.0]
      *               [ 1.0   1.0]
      *
      */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

    /**
     *  Build vector [ 2.0]
     *               [ 4.0]
     */
    double v_data[2] = {2.0, 4.0};
    CAROM::Vector v(v_data, 2, false, true);

    /**
      *  [ 1.0   0.0] [ 2.0] = [ 2.0]
      *  [ 1.0   1.0] [ 4.0]   [ 6.0]
      *
      */
    CAROM::Vector *w = NULL;
    asymmetric_matrix.mult(v, w);
    EXPECT_FALSE(w->distributed());
    EXPECT_EQ(w->dim(), 2);
    EXPECT_DOUBLE_EQ((*w)(0), 2.0);
    EXPECT_DOUBLE_EQ((*w)(1), 6.0);
    delete w;
}

TEST(MatrixSerialTest, Test_mult_Vector_Vector_reference)
{
    /**
      *  Build matrix [ 1.0   0.0]
      *               [ 1.0   1.0]
      *
      */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

    /**
     *  Build vector [ 2.0]
     *               [ 4.0]
     */
    double v_data[2] = {2.0, 4.0};
    CAROM::Vector v(v_data, 2, false, true);

    /**
      *  [ 1.0   0.0] [ 2.0] = [ 2.0]
      *  [ 1.0   1.0] [ 4.0]   [ 6.0]
      *
      */
    CAROM::Vector w(2, false);
    asymmetric_matrix.mult(v, w);
    EXPECT_FALSE(w.distributed());
    EXPECT_EQ(w.dim(), 2);
    EXPECT_DOUBLE_EQ(w(0), 2.0);
    EXPECT_DOUBLE_EQ(w(1), 6.0);
}

TEST(MatrixSerialTest, Test_multPlus)
{
    /**
      *  Build matrix [ 1.0   0.0]
      *               [ 1.0   1.0]
      *
      */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

    /**
     *  Build vector [ 2.0]
     *               [ 4.0]
     */
    double v_data[2] = {2.0, 4.0};
    CAROM::Vector v(v_data, 2, false, true);

    /**
     *  Build vector [ 0.5]
     *               [ 1.0]
     */
    double w_data[2] = {0.5, 1.0};
    CAROM::Vector w(w_data, 2, false, true);

    /**
    *  2 * [ 1.0   0.0] [ 0.5] + [ 2.0] = [ 3.0]
    *      [ 1.0   1.0] [ 1.0]   [ 4.0]   [ 7.0]
    *
    */
    asymmetric_matrix.multPlus(v, w, 2.0);
    EXPECT_DOUBLE_EQ(v(0), 3.0);
    EXPECT_DOUBLE_EQ(v(1), 7.0);
}


TEST(MatrixSerialTest, Test_transposeMult_Vector_reference)
{
    /**
      *  Build matrix [ 1.0   0.0]
      *               [ 1.0   1.0]
      *
      */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

    /**
     *  Build vector [ 2.0]
     *               [ 4.0]
     */
    double v_data[2] = {2.0, 4.0};
    CAROM::Vector v(v_data, 2, false, true);

    /**
      *  [ 1.0   0.0]^{T} [ 2.0] = [ 6.0]
      *  [ 1.0   1.0]     [ 4.0]   [ 4.0]
      *
      */
    CAROM::Vector *w;
    w = asymmetric_matrix.transposeMult(v);
    EXPECT_FALSE(w->distributed());
    EXPECT_EQ(w->dim(), 2);
    EXPECT_DOUBLE_EQ((*w)(0), 6.0);
    EXPECT_DOUBLE_EQ((*w)(1), 4.0);
    delete w;
}

TEST(MatrixSerialTest, Test_transposeMult_Vector_pointer)
{
    /**
      *  Build matrix [ 1.0   0.0]
      *               [ 1.0   1.0]
      *
      */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

    /**
     *  Build vector [ 2.0]
     *               [ 4.0]
     */
    double v_data[2] = {2.0, 4.0};
    CAROM::Vector *v = new CAROM::Vector(v_data, 2, false, true);

    /**
     *  [ 1.0   0.0]^{T} [ 2.0] = [ 6.0]
     *  [ 1.0   1.0]     [ 4.0]   [ 4.0]
     *
     */
    CAROM::Vector *w;
    w = asymmetric_matrix.transposeMult(v);
    EXPECT_FALSE(w->distributed());
    EXPECT_EQ(w->dim(), 2);
    EXPECT_DOUBLE_EQ((*w)(0), 6.0);
    EXPECT_DOUBLE_EQ((*w)(1), 4.0);
    delete w;
    delete v;
}

TEST(MatrixSerialTest, Test_transposeMult_Vector_Vector_pointer)
{
    /**
      *  Build matrix [ 1.0   0.0]
      *               [ 1.0   1.0]
      *
      */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

    /**
     *  Build vector [ 2.0]
     *               [ 4.0]
     */
    double v_data[2] = {2.0, 4.0};
    CAROM::Vector v(v_data, 2, false, true);

    /**
     *  [ 1.0   0.0]^{T} [ 2.0] = [ 6.0]
     *  [ 1.0   1.0]     [ 4.0]   [ 4.0]
     *
     */
    CAROM::Vector *w = NULL;
    asymmetric_matrix.transposeMult(v, w);
    EXPECT_FALSE(w->distributed());
    EXPECT_EQ(w->dim(), 2);
    EXPECT_DOUBLE_EQ((*w)(0), 6.0);
    EXPECT_DOUBLE_EQ((*w)(1), 4.0);
    delete w;
}

TEST(MatrixSerialTest, Test_transposeMult_Vector_Vector_reference)
{
    /**
      *  Build matrix [ 1.0   0.0]
      *               [ 1.0   1.0]
      *
      */
    double asymmetric[4] = {1.0, 0.0, 1.0, 1.0};
    CAROM::Matrix asymmetric_matrix(asymmetric, 2, 2, false, true);

    /**
     *  Build vector [ 2.0]
     *               [ 4.0]
     */
    double v_data[2] = {2.0, 4.0};
    CAROM::Vector v(v_data, 2, false, true);

    /**
     *  [ 1.0   0.0]^{T} [ 2.0] = [ 6.0]
     *  [ 1.0   1.0]     [ 4.0]   [ 4.0]
     *
     */
    CAROM::Vector w(2, false);
    asymmetric_matrix.transposeMult(v, w);
    EXPECT_FALSE(w.distributed());
    EXPECT_EQ(w.dim(), 2);
    EXPECT_DOUBLE_EQ(w(0), 6.0);
    EXPECT_DOUBLE_EQ(w(1), 4.0);
}

TEST(MatrixOuterProduct, Test_outerProduct_serial)
{
    /**
     *  Build vector [ 1.0 ] and vector [ 3.0 ]
     *               [ 2.0 ]            [ 4.0 ]
     *                                  [ 5.0 ]
     */

    double v_data[2] = { 1.0, 2.0 };
    double w_data[3] = { 3.0, 4.0, 5.0 };

    CAROM::Vector v(v_data, 2, false, true);
    CAROM::Vector w(w_data, 3, false, true);

    CAROM::Matrix outer_product_vw = CAROM::outerProduct(v, w);
    EXPECT_EQ(outer_product_vw.numRows(), 2);
    EXPECT_EQ(outer_product_vw.numColumns(), 3);
    EXPECT_FALSE(outer_product_vw.distributed());
    EXPECT_DOUBLE_EQ(outer_product_vw(0, 0), 3.0);
    EXPECT_DOUBLE_EQ(outer_product_vw(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(outer_product_vw(0, 2), 5.0);
    EXPECT_DOUBLE_EQ(outer_product_vw(1, 0), 6.0);
    EXPECT_DOUBLE_EQ(outer_product_vw(1, 1), 8.0);
    EXPECT_DOUBLE_EQ(outer_product_vw(1, 2), 10.0);

    CAROM::Matrix outer_product_wv = CAROM::outerProduct(w, v);
    EXPECT_EQ(outer_product_wv.numRows(), 3);
    EXPECT_EQ(outer_product_wv.numColumns(), 2);
    EXPECT_FALSE(outer_product_wv.distributed());
    EXPECT_DOUBLE_EQ(outer_product_wv(0, 0), 3.0);
    EXPECT_DOUBLE_EQ(outer_product_wv(0, 1), 6.0);
    EXPECT_DOUBLE_EQ(outer_product_wv(1, 0), 4.0);
    EXPECT_DOUBLE_EQ(outer_product_wv(1, 1), 8.0);
    EXPECT_DOUBLE_EQ(outer_product_wv(2, 0), 5.0);
    EXPECT_DOUBLE_EQ(outer_product_wv(2, 1), 10.0);
}

TEST(DiagonalMatrixFactorySerialTest, Test_123vector)
{
    /** Set up the vector [1, 2, 3]^{T} */
    CAROM::Vector w(3, false);
    w(0) = 1.0;
    w(1) = 2.0;
    w(2) = 3.0;

    /**
     *  \diag([1, 2, 3])^{T} = [ 1.0  0.0  0.0 ]
     *                         [ 0.0  2.0  0.0 ]
     *                         [ 0.0  0.0  3.0 ]
     */
    CAROM::Matrix diagonalMatrix = DiagonalMatrixFactory(w);
    EXPECT_TRUE(diagonalMatrix.distributed() == w.distributed());
    EXPECT_TRUE(diagonalMatrix.numRows() == w.dim());
    EXPECT_TRUE(diagonalMatrix.numColumns() == w.dim());
    EXPECT_TRUE(diagonalMatrix.numRows() == diagonalMatrix.numColumns());
    EXPECT_DOUBLE_EQ(diagonalMatrix(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(diagonalMatrix(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(diagonalMatrix(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(diagonalMatrix(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(diagonalMatrix(1, 1), 2.0);
    EXPECT_DOUBLE_EQ(diagonalMatrix(1, 2), 0.0);
    EXPECT_DOUBLE_EQ(diagonalMatrix(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(diagonalMatrix(2, 1), 0.0);
    EXPECT_DOUBLE_EQ(diagonalMatrix(2, 2), 3.0);
}

TEST(IdentityMatrixFactorySerialTest, Test_3vector)
{
    /** Set up the vector [1, 2, 3]^{T} */
    CAROM::Vector w(3, false);
    w(0) = 1.0;
    w(1) = 2.0;
    w(2) = 3.0;

    CAROM::Matrix identityMatrix = IdentityMatrixFactory(w);
    EXPECT_TRUE(identityMatrix.distributed() == w.distributed());
    EXPECT_TRUE(identityMatrix.numRows() == w.dim());
    EXPECT_TRUE(identityMatrix.numColumns() == w.dim());
    EXPECT_TRUE(identityMatrix.numRows() == identityMatrix.numColumns());
    EXPECT_DOUBLE_EQ(identityMatrix(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(identityMatrix(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(identityMatrix(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(identityMatrix(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(identityMatrix(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(identityMatrix(1, 2), 0.0);
    EXPECT_DOUBLE_EQ(identityMatrix(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(identityMatrix(2, 1), 0.0);
    EXPECT_DOUBLE_EQ(identityMatrix(2, 2), 1.0);
}

TEST(MatrixParallelTest, Test_distribute_and_gather)
{
    int is_mpi_initialized, is_mpi_finalized;
    MPI_Initialized(&is_mpi_initialized);
    MPI_Finalized(&is_mpi_finalized);
    if (!is_mpi_initialized) return;

    const MPI_Comm my_comm = MPI_COMM_WORLD;
    int my_rank = -1, num_procs = -1;
    MPI_Comm_size(my_comm, &num_procs);
    MPI_Comm_rank(my_comm, &my_rank);

    int total_rows = 5;
    CAROM::Matrix answer(total_rows, total_rows, false);
    EXPECT_FALSE(answer.distributed());
    for (int i = 0; i < total_rows; i++)
        for (int j = 0; j < total_rows; j++)
            answer.item(i, j) = static_cast<double> (i * j);

    int local_rows = CAROM::split_dimension(total_rows, MPI_COMM_WORLD);
    std::vector<int> row_offsets;
    int total_rows_check = CAROM::get_global_offsets(local_rows, row_offsets, MPI_COMM_WORLD);
    EXPECT_EQ(total_rows, total_rows_check);

    CAROM::Matrix test(answer);
    test.distribute(local_rows);
    for (int local_i = 0, global_i = row_offsets[my_rank]; local_i < local_rows; local_i++, global_i++)
        for (int j = 0; j < answer.numColumns(); j++)
            EXPECT_DOUBLE_EQ(test.item(local_i, j), answer.item(global_i, j));

    test.gather();
    for (int i = 0; i < total_rows; i++)
        for (int j = 0; j < total_rows; j++)
            EXPECT_DOUBLE_EQ(test.item(i, j), answer.item(i, j));

}


int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
#else // #ifndef CAROM_HAS_GTEST
int main()
{
    std::cout << "libROM was compiled without Google Test support, so unit "
              << "tests have been disabled. To enable unit tests, compile "
              << "libROM with Google Test support." << std::endl;
}
#endif // #endif CAROM_HAS_GTEST
