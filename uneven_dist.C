/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 * Written by William Arrighi wjarrighi@llnl.gov
 * CODE-686965
 * All rights reserved.
 *
 * This file is part of libROM.
 * For details, see https://computation.llnl.gov/librom
 * Please also read README_BSD_NOTICE.
 *
 * Redistribution and use in source and binary forms, with or without
 * modifications, are permitted provided that the following conditions are met:
 *
 *    o Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the disclaimer below.
 *    o Redistribution in binary form must reproduce the above copyright
 *      notice, this list of conditions and the disclaimer (as noted below) in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *    o Neither the name of the LLNS/LLNL nor the names of its contributors may
 *      be used to endorse or promote products derived from this software
 *      without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
 * LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OR SUCH DAMAGE.
 *
 *****************************************************************************/

// Description: A test of both the incremental fast update algorithm and
//              sampler in which the complete system in unevenly distributed
//              across the processors.

#include "IncrementalSVDBasisGenerator.h"

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
   // dimension of the problem and the offset into the array of values defining
   // the entire system for this processor.
   int dim;
   int offset;
   if (size == 1) {
      dim = 6;
      offset = 0;
   }
   else if (size == 2) {
      if (rank == 0) {
         dim = 4;
         offset = 0;
      }
      else {
         dim = 2;
         offset = 4;
      }
   }
   else if (size == 3) {
      if (rank == 0) {
         dim = 1;
         offset = 0;
      }
      else if (rank == 1) {
         dim = 2;
         offset = 1;
      }
      else {
         dim = 3;
         offset = 3;
      }
   }
   else if (size == 4) {
      if (rank == 0) {
         dim = 1;
         offset = 0;
      }
      else if (rank == 1) {
         dim = 2;
         offset = 1;
      }
      else if (rank == 2) {
         dim = 2;
         offset = 3;
      }
      else {
         dim = 1;
         offset = 5;
      }
   }
   else if (size == 5) {
      if (rank == 0) {
         dim = 1;
         offset = 0;
      }
      else if (rank == 1) {
         dim = 1;
         offset = 1;
      }
      else if (rank == 2) {
         dim = 1;
         offset = 2;
      }
      else if (rank == 3) {
         dim = 1;
         offset = 3;
      }
      else {
         dim = 2;
         offset = 4;
      }
   }
   else if (size == 6) {
      dim = 1;
      offset = rank;
   }
   else {
      if (rank == 0) {
         printf("Illegal number of procs.\n");
         printf("Allowed number of procs is 1, 2, 3, 4, 5, or 6.\n");
      }
      return 1;
   }

   // Construct the incremental basis generator given the dimension just
   // computed to use the fast update incremental algorithm and the incremental
   // sampler.
   CAROM::IncrementalSVDBasisGenerator inc_basis_generator(dim,
      1.0e-2,
      false,
      true,
      2,
      1.0e-6,
      2,
      1.0e-2,
      0.11,
      "",
      false,
      false,
      false,
      CAROM::Database::HDF5,
      0.1,
      0.8,
      5.0,
      true);

   // Define the values for the first sample.
   double vals0[6] = {1.0, 6.0, 3.0, 8.0, 17.0, 9.0};

   // Define the values for the second sample.
   double vals1[6] = {2.0, 7.0, 4.0, 9.0, 18.0, 10.0};

   bool status = false;

   // Take the first sample.
   if (inc_basis_generator.isNextSample(0.0)) {
      status = inc_basis_generator.takeSample(&vals0[offset], 0.0, 0.11);
      if (status) {
         inc_basis_generator.computeNextSampleTime(&vals0[offset],
            &vals0[offset],
            0.0);
      }
   }

   // Take the second sample.
   if (status && inc_basis_generator.isNextSample(0.11)) {
      status = inc_basis_generator.takeSample(&vals1[offset], 0.11, 0.11);
      if (status) {
         inc_basis_generator.computeNextSampleTime(&vals1[offset],
            &vals1[offset],
            0.11);
      }
   }

   // Finalize MPI and return.
   MPI_Finalize();
   return !status;
}
