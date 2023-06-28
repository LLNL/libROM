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

#include "mpi_utils.h"
#include "Utilities.h"

#include <iomanip>
#include <stdlib.h>
#include <sys/stat.h>

namespace CAROM {

int
split_dimension(const int &dim, const MPI_Comm &comm)
{
    int mpi_init;
    MPI_Initialized(&mpi_init);
    CAROM_VERIFY(mpi_init != 0);

    int d_num_procs, d_rank;
    MPI_Comm_rank(comm, &d_rank);
    MPI_Comm_size(comm, &d_num_procs);

    int local_dim = dim / d_num_procs;
    if (dim % d_num_procs > d_rank)
        local_dim++;

    return local_dim;
}

int
get_global_offsets(const int &local_dim, std::vector<int> &offsets, const MPI_Comm &comm)
{
    int mpi_init;
    MPI_Initialized(&mpi_init);
    CAROM_VERIFY(mpi_init != 0);

    int d_num_procs, d_rank;
    MPI_Comm_rank(comm, &d_rank);
    MPI_Comm_size(comm, &d_num_procs);

    int dim = local_dim;
    CAROM_VERIFY(MPI_Allreduce(MPI_IN_PLACE,
                                &dim,
                                1,
                                MPI_INT,
                                MPI_SUM,
                                comm) == MPI_SUCCESS);

    offsets.resize(d_num_procs + 1);
    offsets[d_num_procs] = dim;
    offsets[d_rank] = local_dim;
    CAROM_VERIFY(MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT,
                               &offsets[0], 1, MPI_INT,
                               comm) == MPI_SUCCESS);
    for (int i = d_num_procs - 1; i >= 0; i--)
        offsets[i] = offsets[i + 1] - offsets[i];

    CAROM_VERIFY(offsets[0] == 0);

    return dim;
}

}
