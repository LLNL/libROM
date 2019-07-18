/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

#include "test-common.hpp"
#include "catch.hpp"

#include "../src/SVD.hpp"
#include "../src/BLAS.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#define NOISY_TESTS false

TEST_CASE("Computing SVD of a preconstructed matrix", "[SVD]")
{
    if (NOISY_TESTS && mpi_rank() == 0)
    {
        std::cout << "Testing SVD with a preconstructed matrix\n";
    }
    // Small blocks so that U is distributed
    constexpr int mb = 2, nb = 2;

    const int totaldim = mpi_size() * 12;
    const double entry = 1 / std::sqrt(3*mpi_size());

    // Create a preconstructed unitary U
    double* Udata = new double[48];
    std::fill(Udata, Udata + 48, 0.0);
    for (int j = 0; j < 4; ++j) {
        Udata[12*j + j] = entry;
        Udata[4+12*j+j] = entry;
        Udata[8+12*j+j] = entry;
    }

    // Create sigmas as 10 ^ (normal random variable)
    std::vector<double> sigmas(4);
    randn(sigmas.data(), 4);
    for (auto& s: sigmas) { s = std::pow(10, s); }
    std::sort(sigmas.begin(), sigmas.end(), std::greater<double>());
    MPI_Bcast(sigmas.data(), 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (NOISY_TESTS && mpi_rank() == 0) {
        std::cout << "Generated sigmas: [";
        for (auto s: sigmas)
            std::cout << s << ' ';
        std::cout << "]\n";
    }

    // Scale columns of U by sigmas.
    for (int j = 0; j < 4; ++j) {
        dscal_wrapper(12, sigmas[j], Udata+j*12, 1);
    }

    // Construct global U.
    ScalaMat U(totaldim, 4, mb, nb);
    for (int rank = 0; rank < mpi_size(); ++rank) {
        U(rowrange(1+rank*12, (rank+1)*12), colrange(1, 4)) =
            LocalMatrix(Udata, 12, 4, rank, false, COL_MAJOR);
    }

    // Factorize it.
    auto info = svd(U);
    // Store the original data back into U, which has been destroyed.
    for (int rank = 0; rank < mpi_size(); ++rank) {
        U(rowrange(1+rank*12, (rank+1)*12), colrange(1, 4)) =
            LocalMatrix(Udata, 12, 4, rank, false, COL_MAJOR);
    }
    delete[] Udata;

    // Check that singular values match.
    assert(info.S.size() == 4);
    if (NOISY_TESTS && mpi_rank() == 0) {
        std::cout << "Computed sigmas: [";
        for (auto s: info.S)
            std::cout << s << ' ';
        std::cout << "]\n";
    }
    for (int i = 0; i < 4; ++i) {
        info.S[i] -= sigmas[i];
    }
    REQUIRE( dnrm2_wrapper(4, info.S.data(), 1) <= 1e-12 );

    // Multiply the computed U, S, and V together and verify that result matches
    double localsdata[16] = { 0 };
    for (int i = 0; i < 4; ++i) { localsdata[i*4+i] = sigmas[i]; }
    ScalaMat S(4, 4, mb, nb, U.context());
    S(rowrange(1, 4), colrange(1, 4)) = LocalMatrix(localsdata, 4, 4, 0, false, COL_MAJOR);
    
    ScalaMat tmp = S.similar();
    ScalaMat Ucomputed = U.similar();

    gemm(false, 1.0, S, false, *info.Vt, 0.0, tmp);
    gemm(false, 1.0, *info.U, false, tmp, 0.0, Ucomputed);
    axpy(false, -1.0, U, Ucomputed);
    double* diffdata = nullptr;
    if (mpi_rank() == 0) {
        diffdata = new double[totaldim*4];
    }
    LocalMatrix(diffdata, totaldim, 4, 0, false, COL_MAJOR) =
        Ucomputed(rowrange(1, totaldim), colrange(1, 4));
    if (mpi_rank() == 0) {
        localsdata[0] = dnrm2_wrapper(totaldim*4, diffdata, 1);
    }
    MPI_Bcast(localsdata, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    delete[] diffdata;

    REQUIRE( localsdata[0] <= 1e-12 );
}

TEST_CASE("Computing SVD of a random matrix", "[SVD]")
{
    constexpr int m = 300, n = 100;
    constexpr int mb = 32, nb = 32;

    ScalaMat U(m, n, mb, nb);
    ScalaMat V(n, n, mb, nb);
    ScalaMat A(m, n, mb, nb);
    ScalaMat Acopy(m, n, mb, nb);

    make_unitary(U);
    make_unitary(V);

    std::vector<double> sigmas(n);
    randn(sigmas.data(), sigmas.size());
    std::sort(sigmas.begin(), sigmas.end(), std::greater<double>());
    for (auto& s: sigmas) { s = std::pow(10, s); }
    MPI_Bcast(sigmas.data(), sigmas.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Create diagonal matrix of singular values.
    // I initialize this memory with enough to hold an m x n matrix to reuse
    // it later.
    std::unique_ptr<double[]> data = nullptr;
    if (mpi_rank() == 0) {
        data = std::unique_ptr<double[]>(new double[m*n]);
        std::fill(data.get(), data.get() + n*n, 0.0);
    }
    LocalMatrix S(data.get(), n, n, 0, false, COL_MAJOR);
    for (int i = 0; i < n; ++i) {
        if (mpi_rank() == 0) {
            data[i*n + i] = sigmas[i];
        }
    }

    // Create A by multiplying U, S, V together.
    ScalaMat tmp(n, n, mb, nb);
    ScalaMat Sglobal(std::move(S), tmp.context());
    gemm(false, 1.0, Sglobal, true, V, 0.0, tmp);
    gemm(false, 1.0, U, false, tmp, 0.0, A);

    // Do the factorization
    Acopy(rowrange(1, m), colrange(1, n)) = A(rowrange(1, m), colrange(1, n));
    auto info = svd(A);

    // Singular values should match.
    axpby_wrapper(sigmas.size(), -1.0, info.S.data(), 1, 1.0, sigmas.data(), 1);
    REQUIRE(dnrm2_wrapper(sigmas.size(), sigmas.data(), 1) <= 1e-12);

    // I don't check U and V explicitly because the signs can be different and
    // it's too much logic for a test. Instead just multiply everything together
    // and see if it matches.
    if (mpi_rank() == 0) {
        std::fill(data.get(), data.get() + n*n, 0.0);
    }
    for (int i = 0; i < n; ++i) {
        if (mpi_rank() == 0) {
            data[i*n + i] = info.S[i];
        }
    }

    Sglobal = ScalaMat(std::move(S), tmp.context());
    gemm(false, 1.0, Sglobal, false, *info.Vt, 0.0, tmp);
    gemm(false, 1.0, *info.U, false, tmp, 0.0, A);
    axpy(false, -1.0, Acopy, A);
    S = LocalMatrix(data.get(), m, n, 0, false, COL_MAJOR);
    S = A(rowrange(1, m), colrange(1, n));
    double dnorm;
    if (mpi_rank() == 0) {
        dnorm = norm(S);
    }
    MPI_Bcast(&dnorm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    REQUIRE( dnorm <= 1e-10 );
}
