/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

#include "MPI_utils.hpp"

#include <cassert>
#include <iostream>

#include "mpi.h"

namespace ScalaWRAP {

namespace {

void mpi_error_check(int err, const char* name)
{
    if (DO_MPI_ERROR_CHECKS) {
        if (err != MPI_SUCCESS) {
            std::cerr << "MPI error code from " << name << ": " << err << '\n';
            throw MPIExcept();
        }
    }
}

}

std::pair<int, int> mpi_global_params()
{
    static int initialized = 0;
    static int nprocs = -1, rank = -1;

    if (nprocs != -1) {
        assert(rank != -1);
        return std::make_pair(nprocs, rank);
    } else {
        assert(rank == -1);
        if (initialized == 0) {
            int err = MPI_Initialized(&initialized);
            mpi_error_check(err, "MPI_Initialized");
            if (initialized != 0) {
                err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
                mpi_error_check(err, "MPI_Comm_size");
                err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                mpi_error_check(err, "MPI_Comm_rank");
                return std::make_pair(nprocs, rank);
            } else {
                return std::make_pair(1, 0);
            }
        } else {
            assert(nprocs != -1 && rank != -1);
            return std::make_pair(nprocs, rank);
        }
    }
}

int mpi_rank()
{
    return mpi_global_params().second;
}

int mpi_size()
{
    return mpi_global_params().first;
}

} /* namespace ScalaWRAP */