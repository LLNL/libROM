/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Utility functions for MPI distribution.

#ifndef included_mpi_utils_h
#define included_mpi_utils_h

#include "CAROM_config.h"
#include "mpi.h"
#include <vector>

namespace CAROM {

/**
 * @brief Distribute the global size dim into MPI processes as equally as possible.
 *
 * @param[in] dim          Input global size.
 * @param[in] comm         MPI communicator. default value MPI_COMM_WORLD.
 * @param[out] local_dim   (Returned value) Local size assigned to the current MPI process.
 */
int
split_dimension(const int dim, const MPI_Comm &comm=MPI_COMM_WORLD);

/**
 * @brief Save integer offsets for each MPI rank under MPI communicator comm,
 *        where each rank as the size of local_dim.
 *        the sum of local_dim over all MPI processes is returned as the total dimension.
 *
 * @param[in] local_dim     Input local dimension specified for each MPI rank.
 * @param[out] offsets      Resulting global integer offsets split by local_dim.
 * @param[in] comm          MPI communicator. default value MPI_COMM_WORLD.
 * @param[out] dim          (Returned value) Global dimension as the sum of all local_dim.
 */
int
get_global_offsets(const int local_dim, std::vector<int> &offsets,
                   const MPI_Comm &comm=MPI_COMM_WORLD);

}

#endif
