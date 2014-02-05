#include "incremental_svd_rom.h"
#include "static_svd_rom.h"

#include "mpi.h"

#include <stdio.h>

int main(int argc, char* argv[])
{
   if (argc != 4) {
      printf("Usage:  random_test dim num_snapshots num_lin_dep_snapshots\n");
      return 1;
   }
   int dim = atoi(argv[1]);
   int num_snapshots = atoi(argv[2]);
   int num_lin_dep_snapshots = atoi(argv[3]);
   int num_lin_indep_snapshots = num_snapshots - num_lin_dep_snapshots;
   MPI_Init(&argc, &argv);
   CAROM::incremental_svd_rom inc_rom(dim, 1.0e-6, false, num_snapshots, 1);
#ifdef DEBUG_ROMS
   CAROM::static_svd_rom static_rom(dim);
#endif
   int size;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   srand(1);
   for (int i = 0; i < dim*num_snapshots*rank; ++i) {
      double random = rand();
      random = random/RAND_MAX;
   }
   double** M = new double* [num_snapshots];
   for (int i = 0; i < num_snapshots; ++i) {
      M[i] = new double [dim];
   }
   for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < num_snapshots; ++j) {
         double random = rand();
         random = random/RAND_MAX;
         M[j][i] = random;
      }
   }
   for (int i = 0; i < dim*num_snapshots*(size-rank-1); ++i) {
      double random = rand();
      random = random/RAND_MAX;
   }
   for (int i = 0; i < num_lin_dep_snapshots; ++i) {
      int col = num_snapshots - i - 1;
      for (int j = 0; j < dim; ++j) {
         M[col][j] = 0;
      }
      for (int j = 0; j < num_lin_indep_snapshots; ++j) {
         double random = rand();
         random = random/RAND_MAX;
         for (int k = 0; k < dim; ++k) {
            M[col][k] += random*M[j][k];
         }
      }
   }
   double start_inc = MPI_Wtime();
   for (int i = 0; i < num_snapshots; ++i) {
      if (inc_rom.isNextSnapshot(0.01*i)) {
         inc_rom.takeSnapshot(M[i], 0.01*i);
         inc_rom.computeNextSnapshotTime(M[i], M[i], 0.01*i);
      }
   }
   const CAROM::Matrix* inc_model = inc_rom.getModel(0.0);
   double stop_inc = MPI_Wtime();
   double incremental_run_time = stop_inc - start_inc;
   double global_incremental_run_time;
   if (size == 1) {
      global_incremental_run_time = incremental_run_time;
   }
   else {
      MPI_Reduce(&incremental_run_time, &global_incremental_run_time, 1,
                 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   }
   if (rank == 0) {
      printf("incremental run time = %g\n", global_incremental_run_time/size);
   }
#ifndef DEBUG_ROMS
   NULL_USE(inc_model);
#else
   double start_static = MPI_Wtime();
   for (int i = 0; i < num_snapshots; ++i) {
      if (static_rom.isNextSnapshot(0.01*i)) {
         static_rom.takeSnapshot(M[i], 0.01*i);
         static_rom.computeNextSnapshotTime(M[i], M[i], 0.01*i);
      }
   }
   const CAROM::Matrix* static_model = static_rom.getModel(0.0);
   double stop_static = MPI_Wtime();
   double static_run_time = stop_static - start_static;
   double global_static_run_time;
   if (size == 1) {
      global_static_run_time = static_run_time;
   }
   else {
      MPI_Reduce(&static_run_time, &global_static_run_time, 1,
                 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   }
   if (rank == 0) {
      printf("static run time = %g\n", global_static_run_time/size);
   }

   // Compute the product of the tranpose of the static model and the
   // incremental model.  This should be a unitary matrix.
   CAROM::Matrix* test = static_model->TransposeMult(*inc_model);
   if (rank == 0) {
      for (int row = 0; row < num_snapshots; ++row) {
         for (int col = 0; col < num_lin_indep_snapshots; ++col) {
            printf("%.16e ", test->item(row, col));
         }
         printf("\n");
      }
   }
   delete test;
#endif
   for (int i = 0; i < num_snapshots; ++i) {
      delete [] M[i];
   }
   delete [] M;
   MPI_Finalize();
   return 0;
}
