/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

#ifdef CAROM_HAS_GTEST

#include<gtest/gtest.h>
#include "linalg/scalapack_wrapper.h"
#include "linalg/BasisGenerator.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring> // for memcpy
#include <random>
#include "mpi.h"
#include "utils/mpi_utils.h"

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
    SUCCEED();
}

TEST(StaticSVDTest, Test_StaticSVD)
{
    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int totaldim = nprocs * 12;
    const double entry = 1 / sqrt(static_cast<double>(3 * nprocs));

    // Columns are constructed to be orthogonal with unit norm.
    std::vector<std::vector<double>> columns;
    for (int i = 0; i < 4; ++i) {
        columns.emplace_back(12u);
        for (int j = i; j < 12; j += 4)
            columns[i][j] = entry;
    }

    std::vector<double> sigmas(4);

    if (rank == 0) {
        std::normal_distribution<double> dist;
        auto generator = std::default_random_engine(std::random_device()());
        for (auto& x: sigmas) {
            x = dist(generator);
            x = x < 0 ? -x : x;
        }
        std::sort(sigmas.begin(), sigmas.end(), std::greater<double>());
    }
    MPI_Bcast(sigmas.data(), 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Generated singular values: ");
        for (int i = 0; i < 4; ++i)
            printf("%8.4E  ", sigmas[i]);
        printf("\n");
    }

    {   /* Wrap the sampler in its own scope; ScaLAPACK will call MPI_Comm_finalize
           in the destructor after MPI_Finalize otherwise. */

        // Construct an SVDSampler to send our matrix. I take V = I for simplicity,
        // so the matrix A that we factor is just the columns scaled by the sigmas.
        CAROM::Options static_svd_options = CAROM::Options(12, 4).
                                            setMaxBasisDimension(4).setDebugMode(true);
        CAROM::BasisGenerator sampler(static_svd_options, false);

        for (unsigned j = 0; j < 4; ++j) {
            std::vector<double> similar(columns[j]);
            for (unsigned i = 0; i < 12; ++i)
                similar[i] *= sigmas[j];
            sampler.takeSample(similar.data(), 0, 0);
        }

        auto distU = sampler.getSpatialBasis();
        if (rank == 0) {
            CAROM::Matrix U(totaldim, 4, false);
            memcpy(&U.item(0, 0), &distU->item(0, 0), 48*sizeof(double));
            for (int src = 1; src < nprocs; ++src)
                MPI_Recv(&U.item(12*src, 0), 48, MPI_DOUBLE, src, MPI_ANY_TAG,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (int i = 0; i < totaldim; ++i)
                for (int j = 0; j < 4; ++j)
                    EXPECT_NEAR(abs(columns[j][i%12]), abs(U.item(i, j)), 1.0e-12);

            printf("Actual U                                        Computed U\n");
            printf("========                                        ==========\n");
            for (int i = 0; i < totaldim; ++i) {
                for (int j = 0; j < 4; ++j)
                    printf("%8.4E ", columns[j][i%12]);
                printf("    ");
                for (int j = 0; j < 4; ++j)
                    printf("%8.4E ", U.item(i, j));
                printf("\n");
            }
        } else {
            MPI_Send(&distU->item(0, 0), 48, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        }

    }

    return;
}

TEST(StaticSVDTest, Test_SLPKTranspose)
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
    constexpr int num_total_cols = 5;
    int d_num_rows = CAROM::split_dimension(num_total_rows, MPI_COMM_WORLD);
    std::vector<int> row_offset(d_num_procs + 1);
    const int total_rows = CAROM::get_global_offsets(d_num_rows, row_offset, MPI_COMM_WORLD);
    EXPECT_EQ(total_rows, num_total_rows);

    double* samples = new double[15] {0.5377, -1.3077, -1.3499,
                                      1.8339, -0.4336, 3.0349,
                                      -2.2588, 0.3426, 0.7254,
                                      0.8622, 3.5784, -0.0631,
                                      0.3188, 2.7694, 0.7147};

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

    SLPK_Matrix *d_samples = new SLPK_Matrix;
    std::vector<int> d_dims;
    d_dims.resize(static_cast<unsigned>(d_num_procs));
    MPI_Allgather(&d_num_rows, 1, MPI_INT, d_dims.data(), 1, MPI_INT, MPI_COMM_WORLD);
    int d_total_dim = 0;
    std::vector<int> d_istarts(static_cast<unsigned>(d_num_procs), 0);

    for (unsigned i = 0; i < static_cast<unsigned>(d_num_procs); ++i) {
        d_total_dim += d_dims[i];
        if (i > 0) {
            d_istarts[i] = d_istarts[i-1] + d_dims[i-1];
        }
    }

    /* NOTE: copied from StaticSVD::StaticSVD */
    int d_nprow = d_num_procs;
    int d_npcol = 1;
    int d_blocksize = d_total_dim / d_nprow;
    if (d_total_dim % d_nprow != 0) {
        d_blocksize += 1;
    }
    initialize_matrix(d_samples, d_total_dim, num_total_cols,
                      d_nprow, d_npcol, d_blocksize, d_blocksize);

    int d_num_samples = 0;
    for (int s = 0; s < 5; s++)
    {
        const double *u_in = &samples[3*s];
        for (int rank = 0; rank < d_num_procs; ++rank) {
            scatter_block(d_samples, d_istarts[static_cast<unsigned>(rank)]+1,
                        s+1, u_in, d_dims[static_cast<unsigned>(rank)],
                        1, rank);
        }
    }
    printf("original\n");
    print_debug_info(d_samples);

    int d_blocksize_tr = num_total_cols / d_nprow;
    if (num_total_cols % d_nprow != 0) {
        d_blocksize_tr += 1;
    }

    SLPK_Matrix *transpose = new SLPK_Matrix;
    initialize_matrix(transpose, num_total_cols, d_total_dim,
                      d_nprow, d_npcol, d_blocksize_tr, d_blocksize_tr);

    for (int rank = 0; rank < d_num_procs; ++rank) {
        transpose_submatrix(transpose, 1, d_istarts[static_cast<unsigned>(rank)]+1,
                            d_samples, d_istarts[static_cast<unsigned>(rank)]+1, 1,
                            d_dims[static_cast<unsigned>(rank)], num_total_cols);
    }

    printf("transposed\n");

    print_debug_info(transpose);

    SVDManager *d_factorizer = new SVDManager;
    svd_init(d_factorizer, transpose);
    d_factorizer->dov = 1;
    factorize(d_factorizer);

    print_debug_info(d_factorizer->U);
    print_debug_info(d_factorizer->V);

    int ncolumns = num_total_rows;
    CAROM::Matrix *d_basis = new CAROM::Matrix(d_dims[d_rank], ncolumns, true);
    CAROM::Matrix *d_basis_right = new CAROM::Matrix(ncolumns, num_total_cols, false);

    for (int rank = 0; rank < d_num_procs; ++rank) {
        // gather_transposed_block does the same as gather_block, but transposes
        // it; here, it is used to go from column-major to row-major order.
        gather_transposed_block(&d_basis_right->item(0, 0), d_factorizer->U,
                                1, 1, num_total_cols, ncolumns, rank);
        // V is computed in the transposed order so no reordering necessary.
        gather_block(&d_basis->item(0, 0), d_factorizer->V,
                    1, d_istarts[static_cast<unsigned>(rank)]+1,
                    ncolumns, d_dims[static_cast<unsigned>(rank)], rank);
    }

    delete d_basis, d_basis_right;
    free_matrix_data(d_samples);
    free_matrix_data(transpose);
    delete d_samples, transpose;
}

TEST(StaticSVDTest, Test_StaticSVDClass)
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
    std::vector<int> row_offset(d_num_procs + 1);
    const int total_rows = CAROM::get_global_offsets(d_num_rows, row_offset, MPI_COMM_WORLD);
    EXPECT_EQ(total_rows, num_total_rows);

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

    CAROM::Options svd_options = CAROM::Options(d_num_rows, 3, 1);
    svd_options.setMaxBasisDimension(num_total_rows);
    svd_options.setDebugMode(true);
    svd_options.setRandomizedSVD(false);
    CAROM::BasisGenerator sampler(svd_options, false);
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

TEST(StaticSVDTest, Test_StaticSVDTranspose)
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
    constexpr int num_total_cols = 5;
    int d_num_rows = CAROM::split_dimension(num_total_rows, MPI_COMM_WORLD);
    std::vector<int> row_offset(d_num_procs + 1);
    const int total_rows = CAROM::get_global_offsets(d_num_rows, row_offset, MPI_COMM_WORLD);
    EXPECT_EQ(total_rows, num_total_rows);

    double* sample1 = new double[5] {0.5377, -1.3077, -1.3499};
    double* sample2 = new double[5] {1.8339, -0.4336, 3.0349};
    double* sample3 = new double[5] {-2.2588, 0.3426, 0.7254};
    double* sample4 = new double[5] {0.8622, 3.5784, -0.0631};
    double* sample5 = new double[5] {0.3188, 2.7694, 0.7147};

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

    CAROM::Options svd_options = CAROM::Options(d_num_rows, num_total_cols, 1);
    svd_options.setMaxBasisDimension(num_total_rows);
    svd_options.setDebugMode(true);
    svd_options.setRandomizedSVD(false);
    CAROM::BasisGenerator sampler(svd_options, false);
    sampler.takeSample(&sample1[row_offset[d_rank]], 0, 0);
    sampler.takeSample(&sample2[row_offset[d_rank]], 0, 0);
    sampler.takeSample(&sample3[row_offset[d_rank]], 0, 0);
    sampler.takeSample(&sample4[row_offset[d_rank]], 0, 0);
    sampler.takeSample(&sample5[row_offset[d_rank]], 0, 0);

    const CAROM::Matrix* d_basis = sampler.getSpatialBasis();
    const CAROM::Matrix* d_basis_right = sampler.getTemporalBasis();
    const CAROM::Vector* sv = sampler.getSingularValues();

    EXPECT_EQ(d_basis_right->numRows(), num_total_cols);
    EXPECT_EQ(d_basis_right->numColumns(), 3);
    EXPECT_EQ(d_basis->numRows(), d_num_rows);
    EXPECT_EQ(d_basis->numColumns(), 3);
    EXPECT_EQ(sv->dim(), 3);

    double* d_basis_vals = d_basis->getData();
    double* d_basis_right_vals = d_basis_right->getData();

    for (int i = 0; i < num_total_cols * num_total_rows; i++) {
        EXPECT_NEAR(abs(d_basis_right_vals[i]),
                    abs(basis_right_true_ans[i]), 1e-7);
    }

    for (int i = 0; i < d_num_rows * num_total_rows; i++) {
        EXPECT_NEAR(abs(d_basis_vals[i]),
                    abs(basis_true_ans[row_offset[d_rank] * num_total_rows + i]), 1e-7);
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