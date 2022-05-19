/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: This source file is a test runner that uses the Google Test
// Framework to run unit tests on the CAROM::S_OPT class.

#include <iostream>

#ifdef CAROM_HAS_GTEST
#include<gtest/gtest.h>
#include <mpi.h>
#include "hyperreduction/S_OPT.h"
#include "linalg/Matrix.h"
#include "utils/CSVDatabase.h"
#define _USE_MATH_DEFINES
#include <cmath>

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
    SUCCEED();
}

TEST(S_OPTSerialTest, Test_S_OPT)
{
    std::vector<std::vector<double>> cols;
    std::vector<double> col1;
    std::vector<double> col2;
    std::vector<double> col3;
    std::vector<double> col4;
    std::vector<double> col5;

    CAROM::CSVDatabase database;
    database.getDoubleVector("../tests/s_opt_data/col1.csv", col1);
    database.getDoubleVector("../tests/s_opt_data/col2.csv", col2);
    database.getDoubleVector("../tests/s_opt_data/col3.csv", col3);
    database.getDoubleVector("../tests/s_opt_data/col4.csv", col4);
    database.getDoubleVector("../tests/s_opt_data/col5.csv", col5);

    cols.push_back(col1);
    cols.push_back(col2);
    cols.push_back(col3);
    cols.push_back(col4);
    cols.push_back(col5);

    int num_rows = 1000;
    int num_cols = 5;
    int num_samples = 10;

    double* orthonormal_mat = new double[num_rows * num_cols];

    // Result of S_OPT (f_basis_sampled_inv)
    double* S_OPT_true_ans = new double[50] {
        -2.34972,  0.0965949, 0.606533, 0.638275, -0.386026,
        -2.44988,  0.430449,  1.07449,  1.25031,  -0.647952,
        -2.54602,  0.803596,  1.43966,  1.59905,  -0.715677,
        -3.37144,  3.97716,   1.47431,  1.69748,  -0.171356,
        -3.49946,  4.15045,   1.03944,  1.29441,   0.25575,
        -3.57175,  1.92897,  -3.25661, -2.99323,   2.22553,
        -3.32192, -0.969017, -4.72273, -2.43004,   0.548499,
        -3.38077, -5.25456,  -2.08931,  5.26275,  -2.44579,
        -3.39278, -4.41453,   4.52803, -0.244015,  4.43685,
        -2.01887, -1.11166,   3.18003, -6.02615,  -4.53692
    };

    int index = 0;
    for (int i = 0; i < num_rows; i++)
    {
        for (int j = 0; j < num_cols; j++)
        {
            orthonormal_mat[index] = cols[j][i];
            index++;
        }
    }

    CAROM::Matrix* u = new CAROM::Matrix(orthonormal_mat, num_rows, num_cols, true);

    double* S_OPT_res = NULL;
    std::vector<int> f_sampled_row(num_samples, 0);
    std::vector<int> f_sampled_row_true_ans{61, 92, 113, 257, 281, 410, 466, 545, 638, 716};
    std::vector<int> f_sampled_rows_per_proc(1, 0);
    CAROM::Matrix f_basis_sampled_inv = CAROM::Matrix(num_samples, num_cols, false);
    CAROM::S_OPT(u, num_cols, f_sampled_row, f_sampled_rows_per_proc, f_basis_sampled_inv, 0, 1, num_samples);

    for (int i = 0; i < num_cols; i++) {
        EXPECT_EQ(f_sampled_row[i], f_sampled_row_true_ans[i]);
    }

    // Compare the norm between the S_OPT result and the true S_OPT answer
    double l2_norm_diff = 0.0;
    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < num_cols; j++) {
            l2_norm_diff += pow(abs(S_OPT_true_ans[i * num_cols + j] - f_basis_sampled_inv.item(i, j)), 2);
        }
    }
    l2_norm_diff = sqrt(l2_norm_diff);

    // Allow for some error due to float rounding
    EXPECT_TRUE(l2_norm_diff < 1e-4);
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
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
