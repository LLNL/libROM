/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "mpi.h"

#include "../src/Context.hpp"

int main(int argc, char* argv[])
{
    int err = MPI_Init(&argc, &argv);
    err = Catch::Session().run(argc, argv);
    ScalaWRAP::finalize_scalawrap();
    
    MPI_Finalize();
    return err;
}