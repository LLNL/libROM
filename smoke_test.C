#include "incremental_svd_rom.h"
#include "stdio.h"

int main(int argc, char* argv[])
{
   if (argc != 2) {
      printf("Usage: smoke_test dim\n");
      return 1;
   }
   int dim = atoi(argv[1]);
   CAROM::incremental_svd_rom inc_rom(&argc, &argv, dim, 1.0e-2, false, 1);
   int size;
   MPI_Comm_size(PETSC_COMM_WORLD, &size);
   int rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   if (dim*size != 6) {
      printf("Illegal dim/number of tasks, guess again.\n");
      return 1;
   }
   double next_snapshot_time;
   double vals0[6] = {1.0, 6.0, 3.0, 8.0, 17.0, 9.0};
   if (inc_rom.isNextSnapshot(0.0)) {
      inc_rom.takeSnapshot(&vals0[dim*rank]);
      next_snapshot_time = inc_rom.computeNextSnapshotTime(&vals0[dim*rank],
                                                           &vals0[dim*rank],
                                                           0.0);
   }
   double vals1[6] = {2.0, 7.0, 4.0, 9.0, 18.0, 10.0};
//   double vals1[6] = {1.0, 6.0, 3.0, 8.0, 17.0, 9.0};
   if (inc_rom.isNextSnapshot(0.11)) {
      inc_rom.takeSnapshot(&vals1[dim*rank]);
      next_snapshot_time = inc_rom.computeNextSnapshotTime(&vals1[dim*rank],
                                                           &vals1[dim*rank],
                                                           0.11);
   }
   return 0;
}
