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
#include "linalg/scalapack_wrapper.h"
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

    bool status = false;

    // Define the values for the first sample.
    double vals0[6] = {1.0, 6.0, 3.0, 8.0, 17.0, 9.0};

    // Define the values for the second sample.
    double vals1[6] = {2.0, 7.0, 4.0, 9.0, 18.0, 10.0};

    // Create first static basis generator for snapshot 1
    CAROM::Options static_svd_options = CAROM::Options(dim,
                                        2).setMaxBasisDimension(2);

    // Create an inner scope so destructors are called when out of scope
    if (true) {

        std::unique_ptr<CAROM::BasisGenerator> static_basis_generator;
        static_basis_generator.reset(new CAROM::BasisGenerator(
                                         static_svd_options, false,
                                         "static_smoke1"));

        // Take the first sample.
        static_basis_generator->takeSample(&vals0[dim*rank],0,0.1);
        std::cout << "Writing sample 1" << std::endl;
        static_basis_generator->writeSnapshot();
        static_basis_generator->endSamples();
        std::cout << "Starting generator 2" <<std::endl;

        // FOM would then close and restart ... create 2nd basis generator
        std::unique_ptr<CAROM::BasisGenerator> static_basis_generator2;
        static_basis_generator2.reset(new CAROM::BasisGenerator(
                                          static_svd_options, false,
                                          "static_smoke2"));

        // Take the second sample.
        static_basis_generator2->takeSample(&vals1[dim*rank],0,0.1);
        static_basis_generator2->writeSnapshot(); // "_snapshot" will be added to the base file name
        static_basis_generator2->endSamples();

        // Files are closed in destructor.
    }

    std::cout << "Starting to compute basis from precomputed snapshots" <<
              std::endl;

    // Create basis using 2 already computed snapshots
    std::unique_ptr<CAROM::BasisGenerator> static_basis_generator3;
    static_basis_generator3.reset(new CAROM::BasisGenerator(
                                      static_svd_options, false,
                                      "static_smoke_final"));

    static_basis_generator3->loadSamples("static_smoke1_snapshot","snapshot");
    static_basis_generator3->loadSamples("static_smoke2_snapshot","snapshot");
    std::cout << "Loaded snapshots" << std::endl;
    static_basis_generator3->writeSnapshot();
    int rom_dim = static_basis_generator3->getSpatialBasis()->numColumns();
    std::cout << "U ROM dimension = " << rom_dim << std::endl;

    static_basis_generator3->endSamples();
    static_basis_generator3 = nullptr;


    // Recreate basis using 1 basis generator as a check
    std::unique_ptr<CAROM::BasisGenerator> static_basis_generator4;
    static_basis_generator4.reset(new CAROM::BasisGenerator(
                                      static_svd_options, false,
                                      "static_smoke_check"));

    static_basis_generator4->takeSample(&vals0[dim*rank],0,0.1);
    static_basis_generator4->takeSample(&vals1[dim*rank],0,0.1);

    static_basis_generator4->endSamples();
    static_basis_generator4 = nullptr;


    // Create basis using 2 already computed bases
    static_basis_generator3.reset(new CAROM::BasisGenerator(
                                      static_svd_options, false,
                                      "static_smoke_final_frombasis"));

    static_basis_generator3->loadSamples("static_smoke1","basis");
    static_basis_generator3->loadSamples("static_smoke2","basis");
    std::cout << "Loaded bases" << std::endl;
    static_basis_generator3->writeSnapshot();
    rom_dim = static_basis_generator3->getSpatialBasis()->numColumns();
    std::cout << "U ROM dimension = " << rom_dim << std::endl;

    static_basis_generator3->endSamples();
    static_basis_generator3 = nullptr;

    // Finalize MPI and return.
    MPI_Finalize();
    return !status;
}
