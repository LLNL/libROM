#include "incremental_svd_rom.h"
#include "static_svd_rom.h"
#include <stdio.h>
#include <mpi.h>

// #define DEBUG

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
   CAROM::incremental_svd_rom inc_rom(dim, 1.0e-12, false, 1);
#ifdef DEBUG
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
      if (inc_rom.isNextSnapshot(0.0)) {
         inc_rom.takeSnapshot(M[i]);
         double next_snapshot_time =
            inc_rom.computeNextSnapshotTime(M[i], M[i], 0.0);
      }
   }
   double* inc_model = inc_rom.getModel();
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
#ifdef DEBUG
   double start_static = MPI_Wtime();
   for (int i = 0; i < num_snapshots; ++i) {
      if (static_rom.isNextSnapshot(0.0)) {
         static_rom.takeSnapshot(M[i]);
         double next_snapshot_time =
            static_rom.computeNextSnapshotTime(M[i], M[i], 0.0);
      }
   }
   double* static_model = static_rom.getModel();
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
   for (int S_col = 0; S_col < num_snapshots; ++S_col) {
      for (int I_col = 0; I_col < num_lin_indep_snapshots; ++I_col) {
         double test = 0.0;
         int S_idx = S_col + num_snapshots*dim*rank;
         int I_idx = I_col;
         for (int entry = 0; entry < dim; ++entry) {
            test += static_model[S_idx]*inc_model[I_idx];
            S_idx += num_snapshots;
            I_idx += num_lin_indep_snapshots;
         }
         double global_test;
         if (size == 1) {
            global_test = test;
         }
         else {
            MPI_Reduce(&test, &global_test, 1, MPI_DOUBLE,
                       MPI_SUM, 0, MPI_COMM_WORLD);
         }
         if (rank == 0) {
            printf("%.16e ", global_test);
         }
      }
      if (rank == 0) {
         printf("\n");
      }
   }
#endif
   delete [] inc_model;
   for (int i = 0; i < num_snapshots; ++i) {
      delete [] M[i];
   }
   delete [] M;
   MPI_Finalize();
   return 0;
}
