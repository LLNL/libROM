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
    database.getDoubleVector("../../tests/s_opt_data/col1.csv", col1);
    database.getDoubleVector("../../tests/s_opt_data/col2.csv", col2);
    database.getDoubleVector("../../tests/s_opt_data/col3.csv", col3);
    database.getDoubleVector("../../tests/s_opt_data/col4.csv", col4);
    database.getDoubleVector("../../tests/s_opt_data/col5.csv", col5);

    cols.push_back(col1);
    cols.push_back(col2);
    cols.push_back(col3);
    cols.push_back(col4);
    cols.push_back(col5);

    int num_rows = 1000;
    int num_cols = 5;
    int num_samples = 10;

    double* orthonormal_mat = new double[num_rows * num_cols];

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

    // // Compare the norm between the S_OPT result and the true S_OPT answer
    // double l2_norm_diff = 0.0;
    // for (int i = 0; i < num_cols; i++) {
    //     for (int j = 0; j < num_cols; j++) {
    //         l2_norm_diff += pow(abs(S_OPT_true_ans[i * num_cols + j] - f_basis_sampled_inv(i, j)), 2);
    //     }
    // }
    // l2_norm_diff = sqrt(l2_norm_diff);
    //
    // // Allow for some error due to float rounding
    // EXPECT_TRUE(l2_norm_diff < 1e-5);
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
