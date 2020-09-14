/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: A test of both the static and incremental fast update
//              algorithms and samplers.  Random numbers are generated as the
//              state vectors in such a way that the global state vector is the
//              same for all processor decompositions when
//              dim * number of processors is constant.

#include "IncrementalSVDBasisGenerator.h"
#include "StaticSVDBasisGenerator.h"

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

   // Set the dimension of the problem and number of samples.
   int dim = 10000;
   int num_samples = 10;

   CAROM::IncrementalSVDBasisGeneratorOptions inc_bg_options;
   inc_bg_options.dim = dim;
   inc_bg_options.linearity_tol = 1.0e-6;
   inc_bg_options.fast_update = true;
   inc_bg_options.max_basis_dimension = num_samples;
   inc_bg_options.initial_dt = 1.0e-2;
   inc_bg_options.samples_per_time_interval = num_samples;
   inc_bg_options.sampling_tol = 1.0e-20;
   inc_bg_options.max_time_between_samples = 10.001;

   // Construct the incremental basis generator to use the fast update
   // incremental algorithm and the incremental sampler.
   CAROM::IncrementalSVDBasisGenerator inc_basis_generator(
     inc_bg_options
   );

   // Initialize random number generator.
   srand(1);

   // Allocate an array for each sample.
   double** M = new double* [num_samples];
   for (int i = 0; i < num_samples; ++i) {
      M[i] = new double [dim];
   }

   // Call the random number generator enough times so that this processor
   // generates it's part of the global sample.
   for (int i = 0; i < dim*num_samples*rank; ++i) {
      double random = rand();
      random = random/RAND_MAX;
   }

   // Fill in the samples.
   for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < num_samples; ++j) {
         double random = rand();
         random = random/RAND_MAX;
         M[j][i] = random;
      }
   }

   MPI_Barrier(MPI_COMM_WORLD);
   double t1 = MPI_Wtime();
   // Take the samples.
   bool status = true;
   int samples_taken = 0;
   for (int i = 0; i < num_samples; ++i) {
      if (inc_basis_generator.isNextSample(0.01*i)) {
         status = inc_basis_generator.takeSample(M[i], 0.01*i, 0.01);
         if (!status) {
            break;
         }
         inc_basis_generator.computeNextSampleTime(M[i], M[i], 0.01*i);
         ++samples_taken;
      }
   }
   if (status) {
      inc_basis_generator.endSamples();
   }
   else {
      if (rank == 0) {
         printf("SVD error\n");
      }
   }
   double t2 = MPI_Wtime();
   double run_time = t2 - t1;
   double global_run_time;
   if (size == 1) {
      global_run_time = run_time;
   }
   else {
      MPI_Reduce(&run_time,
         &global_run_time,
         1,
         MPI_DOUBLE,
         MPI_MAX,
         0,
         MPI_COMM_WORLD);
   }
   if (rank == 0) {
      printf("Wall time = %f\n", global_run_time);
   }

   int foo = inc_basis_generator.getSpatialBasis()->numColumns();
   if (rank == 0) {
      printf("Samples taken = %d\n", samples_taken);
      printf("Number of samples = %d\n", foo);
   }

   // Clean up.
   for (int i = 0; i < num_samples; ++i) {
      delete [] M[i];
   }
   delete [] M;

   // Finalize MPI and return.
   MPI_Finalize();
   return !status;
}
