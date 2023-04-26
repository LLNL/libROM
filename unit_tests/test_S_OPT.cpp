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
    // Get the rank of this process, and the number of processors.
    int mpi_init, d_rank, d_num_procs;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

    int num_total_rows = 100;
    int num_cols = 5;
    int num_samples = 10;

    int num_rows = num_total_rows / d_num_procs;
    if (num_total_rows % d_num_procs > d_rank) {
        num_rows++;
    }
    int *row_offset = new int[d_num_procs + 1];
    row_offset[d_num_procs] = num_total_rows;
    row_offset[d_rank] = num_rows;

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

    std::vector<double*> cols;
    for (int i = 0; i < num_cols; i++)
    {
        double* tmp = new double[num_rows];
        cols.push_back(tmp);
    }

    CAROM::CSVDatabase database;
    database.getDoubleArray("../unit_tests/s_opt_data/col1.csv",  cols[0], num_rows,
                            row_offset[d_rank], 1, 1);
    database.getDoubleArray("../unit_tests/s_opt_data/col2.csv",  cols[1], num_rows,
                            row_offset[d_rank], 1, 1);
    database.getDoubleArray("../unit_tests/s_opt_data/col3.csv",  cols[2], num_rows,
                            row_offset[d_rank], 1, 1);
    database.getDoubleArray("../unit_tests/s_opt_data/col4.csv",  cols[3], num_rows,
                            row_offset[d_rank], 1, 1);
    database.getDoubleArray("../unit_tests/s_opt_data/col5.csv",  cols[4], num_rows,
                            row_offset[d_rank], 1, 1);

    double* orthonormal_mat = new double[num_rows * num_cols];

    // Result of S_OPT (f_basis_sampled_inv)
    double* S_OPT_true_ans = new double[50] {
        -0.5918106,  1.141427,   1.142373,   1.477501,  -0.7002456,
        -1.204534,  1.608694,   0.5777631, -0.8686062,  0.9839489,
        -1.21845,  1.569369,   0.4668517, -0.958101,   0.9462786,
        -1.023713,  0.2118963, -1.029118,  -0.3487639, -0.740034,
        -1.016572,  0.1513296, -1.060517,  -0.2849657, -0.774013,
        -1.009973,  0.0911333, -1.088985,  -0.2200843, -0.8032426,
        -1.107097, -0.8899746, -0.7112781,  1.030203,   0.2788396,
        -1.173239, -1.002324,  -0.1872992,  1.061846,   0.969148,
        -1.040098, -1.004302,   0.8125958,  0.2427455,  0.9231714,
        -0.570251, -0.9721371,  1.327513,  -1.113124,  -1.083388
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
    std::vector<int> f_sampled_row_true_ans{0, 19, 21, 48, 49, 50, 72, 79, 90, 97};
    std::vector<int> f_sampled_rows_per_proc(d_num_procs, 0);
    CAROM::Matrix f_basis_sampled_inv = CAROM::Matrix(num_samples, num_cols, false);
    CAROM::S_OPT(u, num_cols, f_sampled_row, f_sampled_rows_per_proc,
                 f_basis_sampled_inv, d_rank, d_num_procs, num_samples);

    int curr_index = 0;
    for (int i = 1; i < f_sampled_rows_per_proc.size(); i++)
    {
        curr_index += f_sampled_rows_per_proc[i - 1];
        for (int j = curr_index; j < curr_index + f_sampled_rows_per_proc[i]; j++)
        {
            f_sampled_row[j] += row_offset[i];
        }
    }

    for (int i = 0; i < num_cols; i++) {
        EXPECT_EQ(f_sampled_row[i], f_sampled_row_true_ans[i]);
    }

    // Compare the norm between the S_OPT result and the true S_OPT answer
    double l2_norm_diff = 0.0;
    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < num_cols; j++) {
            l2_norm_diff += pow(abs(S_OPT_true_ans[i * num_cols + j] -
                                    f_basis_sampled_inv.item(i, j)), 2);
        }
    }
    l2_norm_diff = sqrt(l2_norm_diff);

    // Allow for some error due to float rounding
    EXPECT_TRUE(l2_norm_diff < 1e-4);
}

TEST(S_OPTSerialTest, Test_S_OPT_less_basis_vectors)
{
    // Get the rank of this process, and the number of processors.
    int mpi_init, d_rank, d_num_procs;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

    int num_total_rows = 100;
    int num_cols = 5;
    int num_basis_vectors = 3;
    int num_samples = 5;

    int num_rows = num_total_rows / d_num_procs;
    if (num_total_rows % d_num_procs > d_rank) {
        num_rows++;
    }
    int *row_offset = new int[d_num_procs + 1];
    row_offset[d_num_procs] = num_total_rows;
    row_offset[d_rank] = num_rows;

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

    std::vector<double*> cols;
    for (int i = 0; i < num_cols; i++)
    {
        double* tmp = new double[num_rows];
        cols.push_back(tmp);
    }

    CAROM::CSVDatabase database;
    database.getDoubleArray("../unit_tests/s_opt_data/col1.csv",  cols[0], num_rows,
                            row_offset[d_rank], 1, 1);
    database.getDoubleArray("../unit_tests/s_opt_data/col2.csv",  cols[1], num_rows,
                            row_offset[d_rank], 1, 1);
    database.getDoubleArray("../unit_tests/s_opt_data/col3.csv",  cols[2], num_rows,
                            row_offset[d_rank], 1, 1);
    database.getDoubleArray("../unit_tests/s_opt_data/col4.csv",  cols[3], num_rows,
                            row_offset[d_rank], 1, 1);
    database.getDoubleArray("../unit_tests/s_opt_data/col5.csv",  cols[4], num_rows,
                            row_offset[d_rank], 1, 1);

    double* orthonormal_mat = new double[num_rows * num_cols];

    // Result of S_OPT (f_basis_sampled_inv)
    double* S_OPT_true_ans = new double[15] {
        -1.433786,  2.925154,   2.209402,
        -1.861405,  0.6665424, -1.247851,
        -1.890357,  0.5250462, -1.304832,
        -1.935721,  0.3061151, -1.365598,
        -2.835807, -3.503679,   1.97353
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
    std::vector<int> f_sampled_row_true_ans{0, 44, 46, 49, 90};
    std::vector<int> f_sampled_rows_per_proc(d_num_procs, 0);
    CAROM::Matrix f_basis_sampled_inv = CAROM::Matrix(num_samples,
                                        num_basis_vectors, false);
    CAROM::S_OPT(u, num_basis_vectors, f_sampled_row, f_sampled_rows_per_proc,
                 f_basis_sampled_inv, d_rank, d_num_procs, num_samples);

    int curr_index = 0;
    for (int i = 1; i < f_sampled_rows_per_proc.size(); i++)
    {
        curr_index += f_sampled_rows_per_proc[i - 1];
        for (int j = curr_index; j < curr_index + f_sampled_rows_per_proc[i]; j++)
        {
            f_sampled_row[j] += row_offset[i];
        }
    }

    for (int i = 0; i < num_basis_vectors; i++) {
        EXPECT_EQ(f_sampled_row[i], f_sampled_row_true_ans[i]);
    }

    // Compare the norm between the S_OPT result and the true S_OPT answer
    double l2_norm_diff = 0.0;
    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < num_basis_vectors; j++) {
            l2_norm_diff += pow(abs(S_OPT_true_ans[i * num_basis_vectors + j] -
                                    f_basis_sampled_inv.item(i, j)), 2);
        }
    }
    l2_norm_diff = sqrt(l2_norm_diff);

    // Allow for some error due to float rounding
    EXPECT_TRUE(l2_norm_diff < 1e-4);
}

TEST(S_OPTSerialTest, Test_S_OPT_init_vector)
{
    // Get the rank of this process, and the number of processors.
    int mpi_init, d_rank, d_num_procs;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

    int num_total_rows = 100;
    int num_cols = 5;
    int num_basis_vectors = 3;
    int num_samples = 5;

    int num_rows = num_total_rows / d_num_procs;
    if (num_total_rows % d_num_procs > d_rank) {
        num_rows++;
    }
    int *row_offset = new int[d_num_procs + 1];
    row_offset[d_num_procs] = num_total_rows;
    row_offset[d_rank] = num_rows;

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

    std::vector<double*> cols;
    for (int i = 0; i < num_cols; i++)
    {
        double* tmp = new double[num_rows];
        cols.push_back(tmp);
    }

    CAROM::CSVDatabase database;
    database.getDoubleArray("../unit_tests/s_opt_data/col1.csv",  cols[0], num_rows,
                            row_offset[d_rank], 1, 1);
    database.getDoubleArray("../unit_tests/s_opt_data/col2.csv",  cols[1], num_rows,
                            row_offset[d_rank], 1, 1);
    database.getDoubleArray("../unit_tests/s_opt_data/col3.csv",  cols[2], num_rows,
                            row_offset[d_rank], 1, 1);
    database.getDoubleArray("../unit_tests/s_opt_data/col4.csv",  cols[3], num_rows,
                            row_offset[d_rank], 1, 1);
    database.getDoubleArray("../unit_tests/s_opt_data/col5.csv",  cols[4], num_rows,
                            row_offset[d_rank], 1, 1);

    double* orthonormal_mat = new double[num_rows * num_cols];

    // Result of S_OPT (f_basis_sampled_inv)
    double* S_OPT_true_ans = new double[15] {
        -1.433786,  2.925154,   2.209402,
        -1.861405,  0.6665424, -1.247851,
        -1.890357,  0.5250462, -1.304832,
        -1.935721,  0.3061151, -1.365598,
        -2.835807, -3.503679,   1.97353
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
    std::vector<int> f_sampled_row_true_ans{0, 44, 46, 49, 90};
    std::vector<int> f_sampled_rows_per_proc(d_num_procs, 0);
    std::vector<int> init_samples;

    // Use just the first true sampled element as an initial sample (90)
    if (row_offset[d_rank] <= 90 && row_offset[d_rank + 1] > 90)
    {
        init_samples.push_back(90 - row_offset[d_rank]);
    }

    CAROM::Matrix f_basis_sampled_inv = CAROM::Matrix(num_samples,
                                        num_basis_vectors, false);
    CAROM::S_OPT(u, num_basis_vectors, f_sampled_row, f_sampled_rows_per_proc,
                 f_basis_sampled_inv, d_rank, d_num_procs, num_samples, &init_samples);

    int curr_index = 0;
    for (int i = 1; i < f_sampled_rows_per_proc.size(); i++)
    {
        curr_index += f_sampled_rows_per_proc[i - 1];
        for (int j = curr_index; j < curr_index + f_sampled_rows_per_proc[i]; j++)
        {
            f_sampled_row[j] += row_offset[i];
        }
    }

    for (int i = 0; i < num_basis_vectors; i++) {
        EXPECT_EQ(f_sampled_row[i], f_sampled_row_true_ans[i]);
    }

    // Compare the norm between the S_OPT result and the true S_OPT answer
    double l2_norm_diff = 0.0;
    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < num_basis_vectors; j++) {
            l2_norm_diff += pow(abs(S_OPT_true_ans[i * num_basis_vectors + j] -
                                    f_basis_sampled_inv.item(i, j)), 2);
        }
    }
    l2_norm_diff = sqrt(l2_norm_diff);

    // Allow for some error due to float rounding
    EXPECT_TRUE(l2_norm_diff < 1e-4);

    // Add the second true sampled element as an initial sample (0)
    if (d_rank == 0)
    {
        init_samples.push_back(0);
    }
    CAROM::S_OPT(u, num_basis_vectors, f_sampled_row, f_sampled_rows_per_proc,
                 f_basis_sampled_inv, d_rank, d_num_procs, num_samples);

    curr_index = 0;
    for (int i = 1; i < f_sampled_rows_per_proc.size(); i++)
    {
        curr_index += f_sampled_rows_per_proc[i - 1];
        for (int j = curr_index; j < curr_index + f_sampled_rows_per_proc[i]; j++)
        {
            f_sampled_row[j] += row_offset[i];
        }
    }

    for (int i = 0; i < num_basis_vectors; i++) {
        EXPECT_EQ(f_sampled_row[i], f_sampled_row_true_ans[i]);
    }

    // Compare the norm between the S_OPT result and the true S_OPT answer
    l2_norm_diff = 0.0;
    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < num_basis_vectors; j++) {
            l2_norm_diff += pow(abs(S_OPT_true_ans[i * num_basis_vectors + j] -
                                    f_basis_sampled_inv.item(i, j)), 2);
        }
    }
    l2_norm_diff = sqrt(l2_norm_diff);

    // Allow for some error due to float rounding
    EXPECT_TRUE(l2_norm_diff < 1e-4);
}

TEST(S_OPTSerialTest, Test_S_OPT_QR)
{
    // Get the rank of this process, and the number of processors.
    int mpi_init, d_rank, d_num_procs;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

    int num_total_rows = 100;
    int num_cols = 5;
    int num_basis_vectors = 3;
    int num_samples = 5;

    int num_rows = num_total_rows / d_num_procs;
    if (num_total_rows % d_num_procs > d_rank) {
        num_rows++;
    }
    int *row_offset = new int[d_num_procs + 1];
    row_offset[d_num_procs] = num_total_rows;
    row_offset[d_rank] = num_rows;

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

    std::vector<double*> cols;
    for (int i = 0; i < num_cols; i++)
    {
        double* tmp = new double[num_rows];
        cols.push_back(tmp);
    }

    CAROM::CSVDatabase database;
    database.getDoubleArray("../unit_tests/s_opt_data/col1.csv",  cols[0], num_rows,
                            row_offset[d_rank], 1, 1);
    database.getDoubleArray("../unit_tests/s_opt_data/col2.csv",  cols[1], num_rows,
                            row_offset[d_rank], 1, 1);
    database.getDoubleArray("../unit_tests/s_opt_data/col3.csv",  cols[2], num_rows,
                            row_offset[d_rank], 1, 1);
    database.getDoubleArray("../unit_tests/s_opt_data/col4.csv",  cols[3], num_rows,
                            row_offset[d_rank], 1, 1);
    database.getDoubleArray("../unit_tests/s_opt_data/col5.csv",  cols[4], num_rows,
                            row_offset[d_rank], 1, 1);

    double* orthonormal_mat = new double[num_rows * num_cols];

    // Result of S_OPT (f_basis_sampled_inv)
    double* S_OPT_true_ans = new double[15] {
        -1.433785, -2.925153, -2.209402,
        -1.861404, -0.6665415, 1.24785,
        -1.890357, -0.5250456, 1.304833,
        -1.93572, -0.3061142, 1.365598,
        -2.835807,  3.503679, -1.97353
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
    std::vector<int> f_sampled_row_true_ans{0, 44, 46, 49, 90};
    std::vector<int> f_sampled_rows_per_proc(d_num_procs, 0);
    std::vector<int> init_samples;
    CAROM::Matrix f_basis_sampled_inv = CAROM::Matrix(num_samples,
                                        num_basis_vectors, false);
    CAROM::S_OPT(u, num_basis_vectors, f_sampled_row, f_sampled_rows_per_proc,
                 f_basis_sampled_inv, d_rank, d_num_procs, num_samples, &init_samples, true);

    int curr_index = 0;
    for (int i = 1; i < f_sampled_rows_per_proc.size(); i++)
    {
        curr_index += f_sampled_rows_per_proc[i - 1];
        for (int j = curr_index; j < curr_index + f_sampled_rows_per_proc[i]; j++)
        {
            f_sampled_row[j] += row_offset[i];
        }
    }

    for (int i = 0; i < num_basis_vectors; i++) {
        EXPECT_EQ(f_sampled_row[i], f_sampled_row_true_ans[i]);
    }

    // Compare the norm between the S_OPT result and the true S_OPT answer
    double l2_norm_diff = 0.0;
    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < num_basis_vectors; j++) {
            l2_norm_diff += pow(abs(S_OPT_true_ans[i * num_basis_vectors + j] -
                                    f_basis_sampled_inv.item(i, j)), 2);
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
