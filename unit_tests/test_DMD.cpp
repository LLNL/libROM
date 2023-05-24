

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
// Framework to run unit tests on the CAROM::DMD class.

#include <iostream>

#ifdef CAROM_HAS_GTEST
#include<gtest/gtest.h>
#include <mpi.h>
#include "algo/DMD.h"
#include "linalg/Vector.h"
#define _USE_MATH_DEFINES
#include <cmath>

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
    SUCCEED();
}

TEST(DMDTest, Test_DMD)
{
    // Get the rank of this process, and the number of processors.
    int mpi_init, d_rank, d_num_procs;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

    int num_total_rows = 5;
    int d_num_rows = num_total_rows / d_num_procs;
    if (num_total_rows % d_num_procs > d_rank) {
        d_num_rows++;
    }
    int *row_offset = new int[d_num_procs + 1];
    row_offset[d_num_procs] = num_total_rows;
    row_offset[d_rank] = d_num_rows;

    MPI_Allgather(MPI_IN_PLACE,
                  1,
                  MPI_INT,
                  row_offset,
                  1,
                  MPI_INT,
                  MPI_COMM_WORLD);

    for (int i = d_num_procs - 1; i >= 0; i--) {
        row_offset[i] = row_offset[i + 1] - row_offset[i];
    }

    double* sample1 = new double[5] {0.5377, 1.8339, -2.2588, 0.8622, 0.3188};
    double* sample2 = new double[5] {-1.3077, -0.4336, 0.3426, 3.5784, 2.7694};
    double* sample3 = new double[5] {-1.3499, 3.0349, 0.7254, -0.0631, 0.7147};
    double* prediction_baseline = new double[5] {-0.4344, -0.0974, 0.0542, 1.2544, 0.9610};

    CAROM::DMD dmd(d_num_rows, 1.0);
    dmd.takeSample(&sample1[row_offset[d_rank]], 0.0);
    dmd.takeSample(&sample2[row_offset[d_rank]], 1.0);
    dmd.takeSample(&sample3[row_offset[d_rank]], 2.0);

    dmd.train(2);
    CAROM::Vector* result = dmd.predict(3.0);

    for (int i = 0; i < d_num_rows; i++) {
        EXPECT_NEAR(result->item(i), prediction_baseline[row_offset[d_rank] + i], 1e-3);
    }

    dmd.save("test_DMD");
    CAROM::DMD dmd_load("test_DMD");
    CAROM::Vector* result_load = dmd_load.predict(3.0);

    for (int i = 0; i < d_num_rows; i++) {
        EXPECT_NEAR(result_load->item(i), prediction_baseline[row_offset[d_rank] + i],
                    1e-3);
    }
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
