/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Simple test of the incremental fast update algorithm and
//              incremental sampler.

#include "StaticSVDBasisGenerator.h"
#include "ScalaMat.hpp"
#include "mpi.h"

#include <stdio.h>

using namespace ScalaWRAP;

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
   
    bool status = false;
    
    // Create basis using 2 already computed bases
    std::unique_ptr<CAROM::SVDBasisGenerator> static_basis_generator4;
//20932
    static_basis_generator4.reset(new CAROM::StaticSVDBasisGenerator(21932,
        300,
        "su2_total",
        900));

    std::cout << "Loading snapshots" << std::endl;
    static_basis_generator4->loadSamples("su2_mach039_snapshot","snapshot");
    //static_basis_generator4->loadSamples("su2_mach040_snapshot","snapshot");
    static_basis_generator4->writeSnapshot();
    std::cout << "Computing SVD" << std::endl;
    int rom_dim = static_basis_generator4->getSpatialBasis()->numColumns();
    std::cout << "U ROM Dimension: " << rom_dim << std::endl;
    static_basis_generator4->endSamples();
    static_basis_generator4 = nullptr;

    // Finalize MPI and return.
    MPI_Finalize();
    return !status;
}

