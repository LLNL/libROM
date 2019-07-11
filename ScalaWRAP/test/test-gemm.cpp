/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

#include "../src/BLAS.hpp"
#include "catch.hpp"
#include "test-common.hpp"

#define NOISY_TESTS false

TEST_CASE("GEMM w/ identity matrix", "[BLAS]")
{
    if (NOISY_TESTS && mpi_rank() == 0) {
        std::cout << "Testing GEMM with an identity matrix\n";
    }
    constexpr int size = 32;
    double localdata[size*size] = { 0 };
    double randomdat[size*size] = { 0 };
    for (int i = 0; i < size; ++i) {
        localdata[i*size + i] = 1;
    }

    ScalaMat I(LocalMatrix(localdata, size, size, 0, false, COL_MAJOR));

    randn(randomdat, size*size);
    ScalaMat A(LocalMatrix(randomdat, size, size, 0, false, COL_MAJOR));

    SECTION("Test that identity times a matrix is the matrix")
    {
        double targetdata[size*size] = { 0 };
        ScalaMat target(LocalMatrix(targetdata, size, size, 0, false, COL_MAJOR));
        gemm(false, 1.0, I, false, A, 0.0, target);
        axpy(false, -1.0, A, target);
        LocalMatrix diff(localdata, size, size, 0);
        diff = target(rowrange(1, size), colrange(1, size));
        double diffnorm;
        if (mpi_rank() == 0) {
            diffnorm = norm(diff);
        }
        MPI_Bcast(&diffnorm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        REQUIRE( diffnorm <= 1e-12 );
    }

    SECTION("Test that identity times a matrix transpose is the matrix transposed")
    {
        double targetdata[size*size] = { 0 };
        ScalaMat target(LocalMatrix(targetdata, size, size, 0, false, COL_MAJOR));
        gemm(false, 1.0, I, true, A, 0.0, target);
        axpy(true, -1.0, A, target);
        LocalMatrix diff(localdata, size, size, 0);
        diff = target(rowrange(1, size), colrange(1, size));
        double diffnorm;
        if (mpi_rank() == 0) {
            diffnorm = norm(diff);
        }
        MPI_Bcast(&diffnorm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        REQUIRE( diffnorm <= 1e-12 );
    }

    SECTION("Test that matrix times identity is the matrix")
    {
        double targetdata[size*size] = { 0 };
        ScalaMat target(LocalMatrix(targetdata, size, size, 0, false, COL_MAJOR));
        gemm(false, 1.0, A, false, I, 0.0, target);
        axpy(false, -1.0, A, target);
        LocalMatrix diff(localdata, size, size, 0);
        diff = target(rowrange(1, size), colrange(1, size));
        double diffnorm;
        if (mpi_rank() == 0) {
            diffnorm = norm(diff);
        }
        MPI_Bcast(&diffnorm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        REQUIRE( diffnorm <= 1e-12 );
    }

    SECTION("Test that matrix transpose times identity is the matrix transposed")
    {
        double targetdata[size*size] = { 0 };
        ScalaMat target(LocalMatrix(targetdata, size, size, 0, false, COL_MAJOR));
        gemm(true, 1.0, A, true, I, 0.0, target);
        axpy(true, -1.0, A, target);
        LocalMatrix diff(localdata, size, size, 0);
        diff = target(rowrange(1, size), colrange(1, size));
        double diffnorm;
        if (mpi_rank() == 0) {
            diffnorm = norm(diff);
        }
        MPI_Bcast(&diffnorm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        REQUIRE( diffnorm <= 1e-12 );
    }
}

/*
 * This test case verifies that the identity $ (A B)^T = B^T A^T $ holds.
 */
TEST_CASE("GEMM w/ transpose matrices", "[BLAS]")
{
    // Let's use some decently large matrices so that they're distributed over
    // processes.
    constexpr int m = 4096, n = 1024, k = 2048;

    // Generate random matrices on one process, for simplicity.
    std::unique_ptr<double[]> Adata, Bdata;
    if (mpi_rank() == 0) {
        Adata = std::unique_ptr<double[]>(new double[m*k]);
        Bdata = std::unique_ptr<double[]>(new double[k*n]);
        randn(Adata.get(), m*k);
        randn(Bdata.get(), k*n);
    }
    
    ScalaMat A(m, k);
    A(rowrange(1, m), colrange(1, k)) = LocalMatrix(Adata.get(), m, k, 0, false, COL_MAJOR);
    ScalaMat B(k, n);
    B(rowrange(1, k), colrange(1, n)) = LocalMatrix(Bdata.get(), k, n, 0, false, COL_MAJOR);

    ScalaMat AB(m, n);
    ScalaMat BtAt(n, m);

    gemm(false, 1.0, A, false, B, 0.0, AB);
    gemm(true, 1.0, B, true, A, 0.0, BtAt);
    axpy(true, -1.0, BtAt, AB);

    double* data = nullptr;
    if (mpi_rank() == 0) {
        data = new double[m*n];
    }
    LocalMatrix diff(data, m, n, 0, true, COL_MAJOR);
    diff = AB(rowrange(1, m), colrange(1, n));

    double diffnorm;
    if (mpi_rank() == 0) {
        diffnorm = norm(diff);
        if (NOISY_TESTS) {
            std::cout << "diffnorm = " << diffnorm << '\n';
        }
    }
    MPI_Bcast(&diffnorm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    REQUIRE( diffnorm <= 1e-10 );
}