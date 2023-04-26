/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Simple test of the incremental fast update algorithm and
//              incremental sampler.

#include "linalg/BasisGenerator.h"

#include "mpi.h"

#include <stdio.h>

int
main(
    int argc,
    char* argv[])
{
    // Initialize MPI and get the number of processors and this processor's
    // rank.
    MPI_Init(&argc, &argv);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Given the number of processors and the rank of this processor set the
    // dimension of the problem.
    int dim;
    if (size == 1) {
        dim = 6;
    }
    else if (size == 2) {
        dim = 3;
    }
    else if (size == 3) {
        dim = 2;
    }
    else if (size == 6) {
        dim = 1;
    }
    else {
        if (rank == 0) {
            printf("Illegal number of procs.\n");
            printf("Allowed number of procs is 1, 2, 3, or 6.\n");
        }
        return 1;
    }

    // Construct the incremental basis generator to use the fast update
    // incremental algorithm and the incremental sampler.
    CAROM::Options svd_options = CAROM::Options(dim, 2).setMaxBasisDimension(2)
                                 .setIncrementalSVD(1.0e-2, 1.0e-6, 1.0e-2, 0.11, true).setDebugMode(true);


    CAROM::BasisGenerator inc_basis_generator(
        svd_options, true
    );

    // Define the values for the first sample.
    double vals0[6] = {1.0, 6.0, 3.0, 8.0, 17.0, 9.0};

    // Define the values for the second sample.
    double vals1[6] = {2.0, 7.0, 4.0, 9.0, 18.0, 10.0};

    bool status = false;

    // Take the first sample.
    if (inc_basis_generator.isNextSample(0.0)) {
        status = inc_basis_generator.takeSample(&vals0[dim*rank], 0.0, 0.11);
        if (status) {
            inc_basis_generator.computeNextSampleTime(&vals0[dim*rank],
                    &vals0[dim*rank],
                    0.0);
        }
    }

    // Take the second sample.
    if (status && inc_basis_generator.isNextSample(0.11)) {
        status = inc_basis_generator.takeSample(&vals1[dim*rank], 0.11, 0.11);
        if (status) {
            inc_basis_generator.computeNextSampleTime(&vals1[dim*rank],
                    &vals1[dim*rank],
                    0.11);
        }
    }

    // Finalize MPI and return.
    MPI_Finalize();
    return !status;
}
