/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: Simple test of the incremental fast update algorithm and
 *              incremental sampler.
 *
 *****************************************************************************/

#include "IncrementalSVDBasisGenerator.h"

#include "mpi.h"

#include <stdio.h>

int main(int argc, char* argv[])
{
   if (argc != 2) {
      printf("Usage: smoke_test dim\n");
      return 1;
   }
   int dim = atoi(argv[1]);
   MPI_Init(&argc, &argv);
   CAROM::IncrementalSVDBasisGenerator inc_basis_generator(dim,
      1.0e-2,
      false,
      2,
      1.0e-2,
      0.11,
      true,
      "");
   int size;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   if (dim*size != 6) {
      printf("Illegal dim/number of tasks, guess again.\n");
      return 1;
   }
   double next_snapshot_time;
   double vals0[6] = {1.0, 6.0, 3.0, 8.0, 17.0, 9.0};
   if (inc_basis_generator.isNextSnapshot(0.0)) {
      inc_basis_generator.takeSnapshot(&vals0[dim*rank], 0.0);
      next_snapshot_time =
         inc_basis_generator.computeNextSnapshotTime(&vals0[dim*rank],
            &vals0[dim*rank],
            0.0);
   }
   double vals1[6] = {2.0, 7.0, 4.0, 9.0, 18.0, 10.0};
   if (inc_basis_generator.isNextSnapshot(0.11)) {
      inc_basis_generator.takeSnapshot(&vals1[dim*rank], 0.11);
      next_snapshot_time =
         inc_basis_generator.computeNextSnapshotTime(&vals1[dim*rank],
            &vals1[dim*rank],
            0.11);
   }
   MPI_Finalize();
   return 0;
}
