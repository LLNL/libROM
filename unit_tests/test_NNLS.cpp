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
#include "linalg/NNLS.h"
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

TEST(NNLS, solve_with_LINF)
{
    int nproc;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int nrow = 50;
    const int ncol = 200;
    const int ncol_local = CAROM::split_dimension(ncol);
    std::vector<int> row_offset(nproc + 1);
    const int total_cols = CAROM::get_global_offsets(ncol_local, row_offset,
                           MPI_COMM_WORLD);
    const double rel_tol = 0.05;
    const double nnls_tol = 1.0e-11;
    const int true_nnz = 44;

    std::default_random_engine generator;
    generator.seed(
        1234); // fix the seed to keep the same result for different nproc.
    std::uniform_real_distribution<> uniform_distribution(0.0, 1.0);
    std::normal_distribution<double> normal_distribution(0.0, 1.0);

    // distribute from a global matrix to keep the same system for different nproc.
    CAROM::Matrix Gt(ncol, nrow, false);
    for (int i = 0; i < ncol; i++)
        for (int j = 0; j < nrow; j++)
            Gt(i, j) = normal_distribution(generator);
    Gt.distribute(ncol_local);

    CAROM::Vector fom_sol(ncol_local, true);
    CAROM::Vector rom_sol(ncol_local, true);
    CAROM::Vector rhs(nrow, false);

    // distribute from a global matrix to keep the same system for different nproc.
    CAROM::Vector fom_sol_serial(ncol, false);
    for (int c = 0; c < ncol; c++)
        fom_sol_serial(c) = uniform_distribution(generator);
    for (int c = 0; c < ncol_local; c++)
        fom_sol(c) = fom_sol_serial(row_offset[rank] + c);

    Gt.transposeMult(fom_sol, rhs);
    rom_sol = 0.0;

    CAROM::Vector rhs_lb(rhs);
    CAROM::Vector rhs_ub(rhs);

    for (int i = 0; i < rhs.dim(); ++i)
    {
        double delta = rel_tol * abs(rhs(i));
        rhs_lb(i) -= delta;
        rhs_ub(i) += delta;
    }

    CAROM::NNLSSolver nnls(nnls_tol, 0, 0, 2);
    nnls.solve_parallel_with_scalapack(Gt, rhs_lb, rhs_ub, rom_sol);

    int nnz = 0;
    for (int i = 0; i < rom_sol.dim(); ++i)
    {
        if (rom_sol(i) != 0.0)
        {
            nnz++;
        }
    }

    std::cout << rank << ": Number of nonzeros in NNLS solution: " << nnz
              << ", out of " << rom_sol.dim() << std::endl;

    MPI_Allreduce(MPI_IN_PLACE, &nnz, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0)
        std::cout << "Global number of nonzeros in NNLS solution: " << nnz << std::endl;
    EXPECT_EQ(nnz, true_nnz);

    // Check residual of NNLS solution
    CAROM::Vector res(Gt.numColumns(), false);
    Gt.transposeMult(rom_sol, res);

    const double normGsol = res.norm();
    const double normRHS = rhs.norm();

    res -= rhs;
    const double relNorm = res.norm() / std::max(normGsol, normRHS);
    std::cout << rank << ": relative residual norm for NNLS solution of Gs = Gw: "
              << relNorm << std::endl;

    double max_error = 0.0;
    for (int k = 0; k < res.dim(); k++)
    {
        max_error = std::max(max_error, abs(res(k) / rhs(k)));
        // printf("rank %d, error(%d): %.3e\n", rank, k, abs(res(k) / rhs(k)));
        EXPECT_TRUE(abs(res(k)) < rel_tol * abs(rhs(k)) + nnls_tol);
    }
    if (rank == 0)
        printf("maximum error: %.5e\n", max_error);
}

TEST(NNLS, solve_with_L2)
{
    int nproc;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int nrow = 30;
    const int ncol = 100;
    const int ncol_local = CAROM::split_dimension(ncol);
    std::vector<int> row_offset(nproc + 1);
    const int total_cols = CAROM::get_global_offsets(ncol_local, row_offset,
                           MPI_COMM_WORLD);
    const double rel_tol = 0.05;
    const double nnls_tol = 1.0e-11;
    const int true_nnz = 21;

    std::default_random_engine generator;
    generator.seed(
        1234); // fix the seed to keep the same result for different nproc.
    std::uniform_real_distribution<> uniform_distribution(0.0, 1.0);
    std::normal_distribution<double> normal_distribution(0.0, 1.0);

    // distribute from a global matrix to keep the same system for different nproc.
    CAROM::Matrix Gt(ncol, nrow, false);
    for (int i = 0; i < ncol; i++)
        for (int j = 0; j < nrow; j++)
            Gt(i, j) = normal_distribution(generator);
    Gt.distribute(ncol_local);

    CAROM::Vector fom_sol(ncol_local, true);
    CAROM::Vector rom_sol(ncol_local, true);
    CAROM::Vector rhs(nrow, false);

    // distribute from a global matrix to keep the same system for different nproc.
    CAROM::Vector fom_sol_serial(ncol, false);
    for (int c = 0; c < ncol; c++)
        fom_sol_serial(c) = uniform_distribution(generator);
    for (int c = 0; c < ncol_local; c++)
        fom_sol(c) = fom_sol_serial(row_offset[rank] + c);

    Gt.transposeMult(fom_sol, rhs);
    rom_sol = 0.0;

    CAROM::Vector rhs_lb(rhs);
    CAROM::Vector rhs_ub(rhs);

    for (int i = 0; i < rhs.dim(); ++i)
    {
        double delta = rel_tol * abs(rhs(i));
        rhs_lb(i) -= delta;
        rhs_ub(i) += delta;
    }

    CAROM::NNLSSolver nnls(nnls_tol, 0, 0, 2, 1.0e-4, 1.0e-14, 100000, 100000,
                           CAROM::NNLS_termination::L2);
    nnls.solve_parallel_with_scalapack(Gt, rhs_lb, rhs_ub, rom_sol);

    int nnz = 0;
    for (int i = 0; i < rom_sol.dim(); ++i)
    {
        if (rom_sol(i) != 0.0)
        {
            nnz++;
        }
    }

    std::cout << rank << ": Number of nonzeros in NNLS solution: " << nnz
              << ", out of " << rom_sol.dim() << std::endl;

    MPI_Allreduce(MPI_IN_PLACE, &nnz, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0)
        std::cout << "Global number of nonzeros in NNLS solution: " << nnz << std::endl;
    EXPECT_EQ(nnz, true_nnz);

    // Check residual of NNLS solution
    CAROM::Vector res(Gt.numColumns(), false);
    Gt.transposeMult(rom_sol, res);

    const double normGsol = res.norm();
    const double normRHS = rhs.norm();

    res -= rhs;
    const double relNorm = res.norm() / std::max(normGsol, normRHS);
    std::cout << rank << ": relative residual norm for NNLS solution of Gs = Gw: "
              << relNorm << std::endl;

    EXPECT_TRUE(res.norm() < rel_tol * normRHS + nnls_tol);
    double max_error = 0.0;
    for (int k = 0; k < res.dim(); k++)
    {
        max_error = std::max(max_error, abs(res(k) / rhs(k)));
    }
    if (rank == 0)
        printf("maximum error: %.5e\n", max_error);
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