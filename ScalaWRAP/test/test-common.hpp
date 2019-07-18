/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

#ifndef TEST_COMMON_HPP
#define TEST_COMMON_HPP

#include "../src/C_interface.h"
#include "../src/LocalMatrix.hpp"

#include <random>

using namespace ScalaWRAP;
using std::cout;

namespace
{

void print_local_mat(const double *dat, int m, int n, MemoryOrder order)
{
    auto index = [m, n, order](int i, int j) {
        if (order == ROW_MAJOR)
        {
            return i * n + j;
        }
        else
        {
            return j * m + i;
        }
    };

    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cout << static_cast<int>(dat[index(i, j)]) << " ";
        }
        cout << '\n';
    }
}

int random_between(int a, int b)
{
    static std::random_device gen;
    return std::uniform_int_distribution<int>(a, b)(gen);
}

std::pair<int, int> get_random_coord(const std::shared_ptr<const Context> &ctxt,
                                     int rank = 0)
{
    int pi = 0, pj = 0;
    if (mpi_rank() == rank)
    {
        int nprow, npcol;
        blacs_gridinfo_wrapper(ctxt->id(), &nprow, &npcol, &pi, &pj);
        pi = random_between(0, nprow-1);
        pj = random_between(0, npcol-1);
    }
    int coords[] = { pi, pj };
    MPI_Bcast(coords, 2, MPI_INT, rank, MPI_COMM_WORLD);
    return std::make_pair(coords[0], coords[1]);
}

std::pair<int, int> best_distribution(int nprocs)
{
    int i = std::floor(std::sqrt(nprocs));
    while (nprocs % i != 0) {
        i -= 1;
    }
    return std::make_pair(i, nprocs / i);
}

} // namespace

namespace ScalaWRAP {

ScalaMat& make_unitary(ScalaMat& U);

void randn(double* dst, int n);

}

#endif /* TEST_COMMON_HPP */