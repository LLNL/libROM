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
    std::vector<double> mat;

    CAROM::CSVDatabase database;
    database.getDoubleVector("s_opt_data/col1.csv", mat, true);
    database.getDoubleVector("s_opt_data/col2.csv", mat, true);
    database.getDoubleVector("s_opt_data/col3.csv", mat, true);
    database.getDoubleVector("s_opt_data/col4.csv", mat, true);
    database.getDoubleVector("s_opt_data/col5.csv", mat, true);

    int num_rows = 1000;
    int num_cols = 5;

    CAROM::Matrix* u = new CAROM::Matrix(mat.data(), num_rows, num_cols, false);
    double* S_OPT_res = NULL;
    // std::vector<int> f_sampled_row(num_cols, 0);
    // std::vector<int> f_sampled_row_true_ans{0, 1, 4, 5, 9};
    // std::vector<int> f_sampled_rows_per_proc(1, 0);
    // CAROM::Matrix f_basis_sampled_inv = CAROM::Matrix(num_cols, num_cols, false);
    CAROM::S_OPT(u, num_cols, true, 0, 1);
    //
    // for (int i = 0; i < num_cols; i++) {
    //     EXPECT_EQ(f_sampled_row[i], f_sampled_row_true_ans[i]);
    // }
    //
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
