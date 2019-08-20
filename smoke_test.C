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
    
   
    bool status = false;
    
    int dim;

    if (rank == 0) { dim = 2616*4;}
    else { dim = 2617*4;}

    // Create basis using 2 already computed bases
    std::unique_ptr<CAROM::SVDBasisGenerator> static_basis_generator4;
    static_basis_generator4.reset(new CAROM::StaticSVDBasisGenerator(dim,
        300,
        "su2_total",
        900));

    std::cout << "Loading snapshots" << std::endl;
    static_basis_generator4->loadSamples("su2_mach039_snapshot","snapshot");
    static_basis_generator4->loadSamples("su2_mach040_snapshot","snapshot");
    //static_basis_generator4->loadSamples("su2_mach039_basis","basis",10);
    //static_basis_generator4->loadSamples("su2_mach040_basis","basis",10);
    std::cout << "Writing snapshots" << std::endl;
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

