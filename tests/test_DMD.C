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
// Framework to run unit tests on the CAROM::DMD class.

#include <iostream>

#ifdef CAROM_HAS_GTEST
#include<gtest/gtest.h>
#include <mpi.h>
#include "../DMD.h"
#include "../Matrix.h"
#define _USE_MATH_DEFINES
#include <cmath>

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
    SUCCEED();
}

TEST(DMDSerialTest, Test_DMD)
{

    // Random input matrix to DMD
    double* rand_mat = new double[25] {
        0.5377,   -1.3077,   -1.3499,   -0.2050,    0.6715,
        1.8339,   -0.4336,    3.0349,   -0.1241,   -1.2075,
       -2.2588,    0.3426,    0.7254,    1.4897,    0.7172,
        0.8622,    3.5784,   -0.0631,    1.4090,    1.6302,
        0.3188,    2.7694,    0.7147,    1.4172,    0.4889
    };

    // // Result of DMD (f_basis_sampled_inv)
    // double* DMD_true_ans = new double[25] {
    //     -0.295811, -0.264874,  1.02179,  -1.05194,  -0.554046,
    //     -0.270643,  1.05349,    0.119162,  0.541832,  0.646459,
    //     -1.33334,  -0.874864,  -0.276067, -0.27327,   0.124747,
    //     0.672776,  0.538704,  -0.735484, -0.794417,  0.388543,
    //     -0.682073, -0.049598, -0.51706,  -0.457748, -1.11295
    // };

    int num_cols = 5;
    int num_rows = 5;

    CAROM::Matrix* u = new CAROM::Matrix(rand_mat, num_rows, num_cols, false);
    CAROM::DMD dmd(u, 0.99, 0, 1);
    // double* DMD_res = NULL;
    // int* f_sampled_row = new int[num_cols] {0};
    // int* f_sampled_row_true_ans = new int[num_cols] {0, 1, 4, 5, 9};
    // int* f_sampled_rows_per_proc = new int[num_cols] {0};
    // CAROM::Matrix f_basis_sampled_inv = CAROM::Matrix(num_cols, num_cols, false);
    // CAROM::DMD(u, num_cols, f_sampled_row, f_sampled_rows_per_proc, f_basis_sampled_inv, 0, 1);
    //
    // for (int i = 0; i < num_cols; i++) {
    //     EXPECT_EQ(f_sampled_row[i], f_sampled_row_true_ans[i]);
    // }
    //
    // // Compare the norm between the DMD result and the true DMD answer
    // double l2_norm_diff = 0.0;
    // for (int i = 0; i < num_cols; i++) {
    //     for (int j = 0; j < num_cols; j++) {
    //         l2_norm_diff += pow(abs(DMD_true_ans[i * num_cols + j] - f_basis_sampled_inv(i, j)), 2);
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
