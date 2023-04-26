/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Simple test of loading precomputed basis or snapshots and
//              computing the static SVD on the loaded samples. Please
//              run in serial. Assumes this file is located in libROM/build.
//              If not, please adjust file address of sample data below.

#include "linalg/BasisGenerator.h"
#include "linalg/scalapack_wrapper.h"
#include "mpi.h"

#include <stdio.h>

int
main(
    int argc,
    char* argv[])
{
    std::string uploaded_data = "snapshot";
    if (argc > 0) {
        if (strcmp(argv[0],"b")
                || strcmp(argv[0],"basis")) std::string uploaded_data = "basis";
    }

    int dim = 6;

    // Create basis using 1 or 2 already computed bases
    std::unique_ptr<CAROM::BasisGenerator> static_basis_generator;

    static_basis_generator.reset(new CAROM::BasisGenerator(
                                     CAROM::Options(dim, 2).setMaxBasisDimension(2), false,
                                     "samples_total"));

    if (uploaded_data == "snapshot") {
        std::cout << "Loading snapshots" << std::endl;
        static_basis_generator->loadSamples("../../tests/load_samples_data/sample1_snapshot",
                                            "snapshot");
        static_basis_generator->loadSamples("../../tests/load_samples_data/sample2_snapshot",
                                            "snapshot");
    }
    else if (uploaded_data == "basis") {
        std::cout << "Loading bases" << std::endl;
        // Load bases. Last input is number of bases to include (allows for truncation)
        static_basis_generator->loadSamples("../../tests/load_samples_data/sample1_basis",
                                            "basis",1);
        static_basis_generator->loadSamples("../../tests/load_samples_data/sample2_basis",
                                            "basis",1);
    }

    std::cout << "Saving data uploaded as a snapshot" << std::endl;
    static_basis_generator->writeSnapshot();

    // Can compute the SVD by calling getSpatialBasis() or endSamples()
    // endSamples() will save the basis file "samples_total..."
    std::cout << "Computing SVD" << std::endl;
    int rom_dim = static_basis_generator->getSpatialBasis()->numColumns();
    std::cout << "U ROM Dimension: " << rom_dim << std::endl;
    static_basis_generator->endSamples();

    static_basis_generator = nullptr;
}
