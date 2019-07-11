/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

/*!
 * @file MPI_utils.hpp
 * 
 * @brief Defines several functions to automate common MPI tasks for the
 * library.
 */

#ifndef MPI_UTILS_HPP
#define MPI_UTILS_HPP

#include <tuple>
#include <utility>

#include "mpi.h"

#define DO_MPI_ERROR_CHECKS false

namespace ScalaWRAP {

struct MPIExcept {};

/*!
 * @brief Retrieve a pair `(nprocs, rank)` corresponding to the common operation
 * of doing `MPI_Comm_size(MPI_COMM_WORLD); MPI_Comm_rank(MPI_COMM_WORLD)`.
 * 
 * Values are cached once computed so that multiple calls do not have to call
 * out to MPI every time.
 */
std::pair<int, int> mpi_global_params();

/*!
 * @brief Get this process's rank in MPI_COMM_WORLD
 */
int mpi_rank();

/*!
 * @brief Get the number of processes in MPI_COMM_WORLD
 */
int mpi_size();

/*!
 * @brief Execute a function on each MPI rank sequentially.
 * 
 * @param[in] f A callable object returning `void` that will be invoked on
 * each process.
 * 
 * @param[in] args... Arguments that are passed to the callable
 */
template <class Func, class... Args>
void mpi_foreach(Func&& f, Args&&... args)
{
    constexpr int foreach_tag = 101;
    int rank, nprocs, dummy = 1;

    MPI_Barrier(MPI_COMM_WORLD);

    std::tie(nprocs, rank) = mpi_global_params();
    const int sender = rank == 0 ? MPI_PROC_NULL : rank - 1;
    const int receiver = rank == nprocs - 1 ? MPI_PROC_NULL : rank + 1;

    MPI_Recv(&dummy, 1, MPI_INT, sender, foreach_tag, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    std::forward<Func>(f)(std::forward<Args>(args)...);
    MPI_Send(&dummy, 1, MPI_INT, receiver, foreach_tag, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
}

} /* namespace ScalaWRAP */

#endif /* MPI_UTILS_HPP */
