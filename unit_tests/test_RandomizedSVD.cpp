

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
// Framework to run unit tests on the CAROM::RandomizedSVD class.

#include <iostream>

#ifdef CAROM_HAS_GTEST
#include<gtest/gtest.h>
#include <mpi.h>
#include "linalg/BasisGenerator.h"
#include "utils/mpi_utils.h"
#define _USE_MATH_DEFINES
#include <cmath>

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
    SUCCEED();
}

TEST(RandomizedSVDTest, Test_RandomizedSVD)
{
    // Get the rank of this process, and the number of processors.
    int mpi_init, d_rank, d_num_procs;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

    constexpr int num_total_rows = 5;
    int d_num_rows = CAROM::split_dimension(num_total_rows, MPI_COMM_WORLD);
    std::vector<int> row_offset;
    CAROM::get_global_offsets(d_num_rows, row_offset, MPI_COMM_WORLD);

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

    CAROM::Options randomized_svd_options = CAROM::Options(d_num_rows, 3, 1);
    randomized_svd_options.setMaxBasisDimension(num_total_rows);
    randomized_svd_options.setDebugMode(true);
    randomized_svd_options.setRandomizedSVD(true);
    CAROM::BasisGenerator sampler(randomized_svd_options, false);
    sampler.takeSample(&sample1[row_offset[d_rank]], 0, 0);
    sampler.takeSample(&sample2[row_offset[d_rank]], 0, 0);
    sampler.takeSample(&sample3[row_offset[d_rank]], 0, 0);

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

TEST(RandomizedSVDTest, Test_RandomizedSVDTransposed)
{

    // Get the rank of this process, and the number of processors.
    int mpi_init, d_rank, d_num_procs;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

    constexpr int num_total_rows = 3;
    constexpr int num_samples = 5;
    int d_num_rows = CAROM::split_dimension(num_total_rows, MPI_COMM_WORLD);
    std::vector<int> row_offset;
    CAROM::get_global_offsets(d_num_rows, row_offset, MPI_COMM_WORLD);

    double* sample1 = new double[3] {0.5377, -1.3077, -1.3499};
    double* sample2 = new double[3] {1.8339, -0.4336, 3.0349};
    double* sample3 = new double[3] {-2.2588, 0.3426, 0.7254};
    double* sample4 = new double[3] {0.8622, 3.5784, -0.0631};
    double* sample5 = new double[3] {0.3188, 2.7694, 0.7147};

    double* basis_right_true_ans = new double[15] {
        3.08158946098238906153E-01,	    -9.49897947980619661301E-02,	-4.50691774108525788911E-01,	
        -1.43697905723455976457E-01,	9.53289043424090820622E-01,	    8.77767692937209131898E-02,	
        -2.23655845793717528158E-02,	-2.10628953513210204207E-01,	8.42235962392685943989E-01,	
        -7.29903965154318323805E-01,	-1.90917141788945754488E-01,	-2.77280930877637610266E-01,	
        -5.92561353877168350834E-01,	-3.74570084880578441089E-02,	5.40928141934190823137E-02,	
    };

    double* basis_true_ans = new double[9] {
        -1.78651649346571794741E-01,	5.44387957786310106023E-01,	    -8.19588518467042281834E-01,	
        -9.49719639253861602768E-01,	-3.13100149275943651084E-01,	-9.50441422536040881122E-04,	
        -2.57130696341890396805E-01,	7.78209514167382598870E-01,	    5.72951792961765460355E-01,
    };

    double* sv_true_ans = new double[3] {
        4.84486375065219387892E+00,     3.66719976398777269821E+00,     2.69114625366671811335E+00,
    };

    CAROM::Options randomized_svd_options = CAROM::Options(d_num_rows, 5, 1);
    randomized_svd_options.setMaxBasisDimension(num_total_rows);
    randomized_svd_options.setDebugMode(true);
    randomized_svd_options.setRandomizedSVD(true);
    CAROM::BasisGenerator sampler(randomized_svd_options, false);
    sampler.takeSample(&sample1[row_offset[d_rank]], 0, 0);
    sampler.takeSample(&sample2[row_offset[d_rank]], 0, 0);
    sampler.takeSample(&sample3[row_offset[d_rank]], 0, 0);
    sampler.takeSample(&sample4[row_offset[d_rank]], 0, 0);
    sampler.takeSample(&sample5[row_offset[d_rank]], 0, 0);

    const CAROM::Matrix* d_basis = sampler.getSpatialBasis();
    const CAROM::Matrix* d_basis_right = sampler.getTemporalBasis();
    const CAROM::Vector* sv = sampler.getSingularValues();

    EXPECT_EQ(d_basis_right->numRows(), num_samples);
    EXPECT_EQ(d_basis_right->numColumns(), 3);
    EXPECT_EQ(d_basis->numRows(), d_num_rows);
    EXPECT_EQ(d_basis->numColumns(), 3);
    EXPECT_EQ(sv->dim(), 3);

    double* d_basis_vals = d_basis->getData();
    double* d_basis_right_vals = d_basis_right->getData();

    for (int i = 0; i < num_samples * 3; i++) {
        EXPECT_NEAR(abs(d_basis_right_vals[i]),
                    abs(basis_right_true_ans[i]), 1e-7);
    }

    for (int i = 0; i < d_num_rows * 3; i++) {
        EXPECT_NEAR(abs(d_basis_vals[i]), abs(basis_true_ans[row_offset[d_rank] * 3 + i]), 1e-7);
    }

    for (int i = 0; i < 3; i++) {
        EXPECT_NEAR(sv->item(i), sv_true_ans[i], 1e-7);
    }
}

TEST(RandomizedSVDTest, Test_RandomizedSVDSmallerSubspace)
{

    // Get the rank of this process, and the number of processors.
    int mpi_init, d_rank, d_num_procs;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

    constexpr int num_total_rows = 5;
    int d_num_rows = CAROM::split_dimension(num_total_rows, MPI_COMM_WORLD);
    std::vector<int> row_offset;
    CAROM::get_global_offsets(d_num_rows, row_offset, MPI_COMM_WORLD);

    double* sample1 = new double[5] {0.5377, 1.8339, -2.2588, 0.8622, 0.3188};
    double* sample2 = new double[5] {-1.3077, -0.4336, 0.3426, 3.5784, 2.7694};
    double* sample3 = new double[5] {-1.3499, 3.0349, 0.7254, -0.0631, 0.7147};

    double* basis_true_ans = new double[10] {
        2.46707456258500545943E-01,	    2.40668223895116051292E-01,	
        1.22600535968995782987E-02,	    6.13241535464043252546E-01,	
        2.43268241912782956504E-02,	    -7.45453669205227487105E-01,	
        -7.73110544854238135315E-01,	9.33140592900994769732E-02,	
        -5.83689483520313912024E-01,	-4.00616848889396720557E-02,	
    };

    double* basis_right_true_ans = new double[4] {
        -1.56565710761715548571E-01,	9.35014519663929455362E-01,	
        -9.78462684973472551775E-01,	-1.90716931883349205545E-01,	
    };

    double* sv_true_ans = new double[2] {
        4.80607940538476441361E+00,	    3.21443716375044896694E+00,
    };

    CAROM::Options randomized_svd_options = CAROM::Options(d_num_rows, 3, 1);
    randomized_svd_options.setMaxBasisDimension(num_total_rows);
    randomized_svd_options.setDebugMode(true);
    randomized_svd_options.setRandomizedSVD(true, 2);
    CAROM::BasisGenerator sampler(randomized_svd_options, false);
    sampler.takeSample(&sample1[row_offset[d_rank]], 0, 0);
    sampler.takeSample(&sample2[row_offset[d_rank]], 0, 0);
    sampler.takeSample(&sample3[row_offset[d_rank]], 0, 0);

    const CAROM::Matrix* d_basis = sampler.getSpatialBasis();
    const CAROM::Matrix* d_basis_right = sampler.getTemporalBasis();
    const CAROM::Vector* sv = sampler.getSingularValues();

    EXPECT_EQ(d_basis->numRows(), d_num_rows);
    EXPECT_EQ(d_basis->numColumns(), 2);
    EXPECT_EQ(d_basis_right->numRows(), 2);
    EXPECT_EQ(d_basis_right->numColumns(), 2);
    EXPECT_EQ(sv->dim(), 2);

    double* d_basis_vals = d_basis->getData();
    double* d_basis_right_vals = d_basis_right->getData();

    for (int i = 0; i < d_num_rows * 2; i++) {
        EXPECT_NEAR(abs(d_basis_vals[i]),
                    abs(basis_true_ans[row_offset[d_rank] * 2 + i]), 1e-7);
    }

    for (int i = 0; i < 4; i++) {
        EXPECT_NEAR(abs(d_basis_right_vals[i]), abs(basis_right_true_ans[i]), 1e-7);
    }

    for (int i = 0; i < 2; i++) {
        EXPECT_NEAR(sv->item(i), sv_true_ans[i], 1e-7);
    }
}

TEST(RandomizedSVDTest, Test_RandomizedSVDTransposedSmallerSubspace)
{

    // Get the rank of this process, and the number of processors.
    int mpi_init, d_rank, d_num_procs;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

    constexpr int num_total_rows = 3;
    constexpr int num_samples = 5;
    constexpr int reduced_rows = 2;
    int d_num_rows = CAROM::split_dimension(num_total_rows, MPI_COMM_WORLD);
    std::vector<int> row_offset;
    CAROM::get_global_offsets(d_num_rows, row_offset, MPI_COMM_WORLD);

    double* sample1 = new double[5] {0.5377, -1.3077, -1.3499};
    double* sample2 = new double[5] {1.8339, -0.4336, 3.0349};
    double* sample3 = new double[5] {-2.2588, 0.3426, 0.7254};
    double* sample4 = new double[5] {0.8622, 3.5784, -0.0631};
    double* sample5 = new double[5] {0.3188, 2.7694, 0.7147};

    double* basis_right_true_ans = new double[10] {
        2.46707456258500545943E-01,	    2.40668223895116051292E-01,	
        1.22600535968995782987E-02,	    6.13241535464043252546E-01,	
        2.43268241912782956504E-02,	    -7.45453669205227487105E-01,	
        -7.73110544854238135315E-01,	9.33140592900994769732E-02,	
        -5.83689483520313912024E-01,	-4.00616848889396720557E-02,	
    };

    double* basis_true_ans = new double[4] {
        -1.56565710761715548571E-01,	9.35014519663929455362E-01,	
        -9.78462684973472551775E-01,	-1.90716931883349205545E-01,	
    };

    double* sv_true_ans = new double[2] {
        4.80607940538476441361E+00,	    3.21443716375044896694E+00,
    };

    CAROM::Options randomized_svd_options = CAROM::Options(d_num_rows, 5, 1);
    randomized_svd_options.setMaxBasisDimension(num_total_rows);
    randomized_svd_options.setDebugMode(true);
    randomized_svd_options.setRandomizedSVD(true, reduced_rows);
    CAROM::BasisGenerator sampler(randomized_svd_options, false);
    sampler.takeSample(&sample1[row_offset[d_rank]], 0, 0);
    sampler.takeSample(&sample2[row_offset[d_rank]], 0, 0);
    sampler.takeSample(&sample3[row_offset[d_rank]], 0, 0);
    sampler.takeSample(&sample4[row_offset[d_rank]], 0, 0);
    sampler.takeSample(&sample5[row_offset[d_rank]], 0, 0);

    const CAROM::Matrix* d_basis = sampler.getSpatialBasis();
    const CAROM::Matrix* d_basis_right = sampler.getTemporalBasis();
    const CAROM::Vector* sv = sampler.getSingularValues();

    d_num_rows = CAROM::split_dimension(reduced_rows, MPI_COMM_WORLD);
    int reduced_rows_check = CAROM::get_global_offsets(d_num_rows, row_offset, MPI_COMM_WORLD);

    EXPECT_EQ(d_basis_right->numRows(), num_samples);
    EXPECT_EQ(d_basis_right->numColumns(), reduced_rows);
    EXPECT_EQ(d_basis->numRows(), d_num_rows);
    EXPECT_EQ(d_basis->numColumns(), reduced_rows);
    EXPECT_EQ(sv->dim(), reduced_rows);

    double* d_basis_vals = d_basis->getData();
    double* d_basis_right_vals = d_basis_right->getData();

    for (int i = 0; i < num_samples * reduced_rows; i++) {
        EXPECT_NEAR(abs(d_basis_right_vals[i]),
                    abs(basis_right_true_ans[i]), 1e-7);
    }

    for (int i = 0; i < d_num_rows * reduced_rows; i++) {
        EXPECT_NEAR(abs(d_basis_vals[i]), abs(basis_true_ans[row_offset[d_rank] * reduced_rows + i]), 1e-7);
    }

    for (int i = 0; i < reduced_rows; i++) {
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
