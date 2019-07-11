/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

#include "test-common.hpp"
#include "../src/C_interface.h"

#include <cassert>
#include <random>

namespace {

std::default_random_engine gen;

void init_gen()
{
    gen.seed(std::random_device()());
}

}

namespace ScalaWRAP {

void randn(double* dst, int n)
{
    init_gen();
    std::normal_distribution<double> dist;
    for (int i = 0; i < n; ++i) {
        dst[i] = dist(gen);
    }
}

ScalaMat& make_unitary(ScalaMat& U)
{
    assert(U.m() >= U.n());
    std::unique_ptr<double[]> data;
    if (mpi_rank() == 0) {
        data = std::unique_ptr<double[]>(new double[U.m() * U.n()]);
    }
    LocalMatrix loc(data.get(), U.m(), U.n(), 0, false, COL_MAJOR);

    if (mpi_rank() == 0) {
        randn(data.get(), U.m()*U.n());
        
        for (int i = 0; i < U.n(); ++i) {
            for (int j = 0; j < i; ++j) {
                double inner = ddot_wrapper(U.m(), data.get() + j*U.m(), 1, 
                                            data.get() + i*U.m(), 1);
                axpby_wrapper(U.m(), -inner, data.get() + j*U.m(), 1, 1.0,
                              data.get() + i*U.m(), 1);
            }
            double norm = dnrm2_wrapper(U.m(), data.get() + i*U.m(), 1);
            dscal_wrapper(U.m(), 1 / norm, data.get() + i*U.m(), 1);
        }
    }

    U(rowrange(1, U.m()), colrange(1, U.n())) = std::move(loc);
    return U;
}

}