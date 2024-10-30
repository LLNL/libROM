/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

#include <iostream>

#ifdef CAROM_HAS_GTEST

#include <gtest/gtest.h>
#include "linalg/scalapack_wrapper.h"
#include "linalg/BasisGenerator.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring> // for memcpy
#include <memory>
#include <random>
#include "mpi.h"
#include "utils/mpi_utils.h"

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
    SUCCEED();
}

TEST(QRfactorizeTest, Test_QR)
{
    // Get the rank of this process, and the number of processors.
    int mpi_init, rank, num_procs;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    constexpr int num_total_rows = 5;
    constexpr int num_columns = 3;
    int loc_num_rows = CAROM::split_dimension(num_total_rows, MPI_COMM_WORLD);
    std::vector<int> row_offset(num_procs + 1);
    const int total_rows = CAROM::get_global_offsets(loc_num_rows, row_offset,
                           MPI_COMM_WORLD);
    EXPECT_EQ(total_rows, num_total_rows);

    double* q_data = new double[15] {
        3.08158946098238906153E-01,     -9.49897947980619661301E-02,     -4.50691774108525788911E-01,
        -1.43697905723455976457E-01,      9.53289043424090820622E-01,      8.77767692937209131898E-02,
        -2.23655845793717528158E-02,     -2.10628953513210204207E-01,      8.42235962392685943989E-01,
        -7.29903965154318323805E-01,     -1.90917141788945754488E-01,     -2.77280930877637610266E-01,
        -5.92561353877168350834E-01,     -3.74570084880578441089E-02,      5.40928141934190823137E-02
    };

    double* r_data = new double[9] {
        -1.78651649346571794741E-01,      5.44387957786310106023E-01,     -8.19588518467042281834E-01,
            0.0,     -3.13100149275943651084E-01,     -9.50441422536040881122E-04,
            0.0,                             0.0,      5.72951792961765460355E-01
        };

    CAROM::Matrix exactQ(loc_num_rows, num_columns, true);
    CAROM::Matrix exactR(r_data, num_columns, num_columns, false);

    for (int i = 0; i < loc_num_rows; i++) {
        for (int j = 0; j < num_columns; j++) {
            exactQ(i, j) = q_data[((i + row_offset[rank]) * num_columns) + j];
        }
    }

    EXPECT_EQ(exactQ.numRows(), loc_num_rows);
    EXPECT_EQ(exactQ.numColumns(), num_columns);

    // Verify that the columns of Q are orthonormal.
    {
        std::unique_ptr<CAROM::Matrix> id = exactQ.transposeMult(exactQ);

        EXPECT_EQ(id->numRows(), num_columns);
        EXPECT_EQ(id->numColumns(), num_columns);

        double maxError = 0.0;
        for (int i = 0; i < num_columns; i++) {
            for (int j = 0; j < num_columns; j++) {
                const double delta_ij = i == j ? 1.0 : 0.0;
                const double error = std::abs((*id)(i,j) - delta_ij);
                maxError = std::max(maxError, error);
            }
        }

        EXPECT_NEAR(maxError, 0.0, 1.0e-15);
    }

    std::unique_ptr<CAROM::Matrix> A = exactQ.mult(exactR);  // Compute A = QR

    // Factorize A
    std::vector<std::unique_ptr<CAROM::Matrix>> QRfactors;
    A->qr_factorize(QRfactors);
    std::unique_ptr<CAROM::Matrix> & Q = QRfactors[0];
    std::unique_ptr<CAROM::Matrix> & R = QRfactors[1];

    EXPECT_EQ(Q->numRows(), loc_num_rows);
    EXPECT_EQ(Q->numColumns(), num_columns);

    // Verify that Q == -exactQ and R == -exactR

    double maxError = 0.0;
    for (int i = 0; i < loc_num_rows; i++) {
        for (int j = 0; j < num_columns; j++) {
            const double error = std::abs(exactQ(i, j) + (*Q)(i, j));
            maxError = std::max(maxError, error);
        }
    }

    EXPECT_NEAR(maxError, 0.0, 1.0e-15);

    EXPECT_EQ(R->numRows(), num_columns);
    EXPECT_EQ(R->numColumns(), num_columns);

    for (int i = 0; i < num_columns; i++) {
        for (int j = 0; j < num_columns; j++) {
            const double error = std::abs(exactR(i, j) + (*R)(i, j));
            maxError = std::max(maxError, error);
        }
    }

    EXPECT_NEAR(maxError, 0.0, 1.0e-15);
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
