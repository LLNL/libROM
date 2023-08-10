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
// Framework to run unit tests on the CAROM::IncrementalSVD class.

#include <iostream>

#ifdef CAROM_HAS_GTEST
#include<gtest/gtest.h>
#include <mpi.h>
#include "linalg/BasisGenerator.h"

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
    SUCCEED();
}

TEST(IncrementalSVDBrandTest, Test_IncrementalSVDBrand)
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

    double* basis_true_ans = new double[15] {
        3.08158946098238906153E-01,      -9.49897947980619661301E-02,      -4.50691774108525788911E-01,
        -1.43697905723455976457E-01,     9.53289043424090820622E-01,      8.77767692937209131898E-02,
        -2.23655845793717528158E-02,     -2.10628953513210204207E-01,     8.42235962392685943989E-01,
        -7.29903965154318323805E-01,     -1.90917141788945754488E-01,     -2.77280930877637610266E-01,
        -5.92561353877168350834E-01,     -3.74570084880578441089E-02,     5.40928141934190823137E-02
    };

    double* basis_right_true_ans = new double[9] {
        -1.78651649346571794741E-01,     5.44387957786310106023E-01,      -8.19588518467042281834E-01,
        -9.49719639253861602768E-01,     -3.13100149275943651084E-01,     -9.50441422536040881122E-04,
        -2.57130696341890396805E-01,     7.78209514167382598870E-01,      5.72951792961765460355E-01
    };

    double* sv_true_ans = new double[3] {
        4.84486375065219387892E+00,      3.66719976398777269821E+00,      2.69114625366671811335E+00
    };

    bool fast_update = true;
    bool fast_update_brand = true;
    CAROM::Options incremental_svd_options = CAROM::Options(d_num_rows, 3, -1, true)
	    .setMaxBasisDimension(num_total_rows)
            .setIncrementalSVD(1e-1,
			       1e-1,
			       1e-1,
			       1e-1,
			       fast_update,
			       fast_update_brand,
			       false);

    CAROM::BasisGenerator sampler(
        incremental_svd_options,
	true,
        "irrelevant.txt");
    sampler.takeSample(&sample1[row_offset[d_rank]], 0, 1e-1);
    sampler.takeSample(&sample2[row_offset[d_rank]], 0, 1e-1);
    sampler.takeSample(&sample3[row_offset[d_rank]], 0, 1e-1);

    const CAROM::Matrix* d_basis = sampler.getSpatialBasis();
    const CAROM::Matrix* d_basis_right = sampler.getTemporalBasis();
    const CAROM::Vector* sv = sampler.getSingularValues();

    EXPECT_EQ(d_basis->numRows(), d_num_rows);
    EXPECT_EQ(d_basis->numColumns(), 3);
    EXPECT_EQ(d_basis_right->numRows(), 3);
    EXPECT_EQ(d_basis_right->numColumns(), 3);
    EXPECT_EQ(sv->dim(), 3);

    double* d_basis_vals = d_basis->getData();
    double* d_basis_right_vals = d_basis_right->getData();


    for (int i = 0; i < d_num_rows * 3; i++) {
        EXPECT_NEAR(abs(d_basis_vals[i]),
                    abs(basis_true_ans[row_offset[d_rank] * 3 + i]), 1e-7);
    }

    for (int i = 0; i < 9; i++) {
        EXPECT_NEAR(abs(d_basis_right_vals[i]), abs(basis_right_true_ans[i]), 1e-7);
    }

    for (int i = 0; i < 3; i++) {
        EXPECT_NEAR(sv->item(i), sv_true_ans[i], 1e-7);
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
