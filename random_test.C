#include "incremental_svd_rom.h"
#include "static_svd_rom.h"
#include "stdio.h"

int main(int argc, char* argv[])
{
   if (argc != 4) {
      printf("Usage:  random_test dim num_snapshots num_lin_dep_snapshots\n");
      return 1;
   }
   int dim = atoi(argv[1]);
   int num_snapshots = atoi(argv[2]);
   int num_lin_dep_snapshots = atoi(argv[3]);
   CAROM::incremental_svd_rom inc_rom(&argc, &argv, dim, 1.0e-12, false, 1);
   CAROM::static_svd_rom static_rom(&argc, &argv, dim);
   int size;
   MPI_Comm_size(PETSC_COMM_WORLD, &size);
   int rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
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
      for (int j = 0; j < num_snapshots - num_lin_dep_snapshots; ++j) {
         double random = rand();
         random = random/RAND_MAX;
         for (int k = 0; k < dim; ++k) {
            M[col][k] += random*M[j][k];
         }
      }
   }
   for (int i = 0; i < num_snapshots; ++i) {
      if (inc_rom.isNextSnapshot(0.0)) {
         inc_rom.takeSnapshot(M[i]);
         double next_snapshot_time =
            inc_rom.computeNextSnapshotTime(M[i], M[i], 0.0);
      }
      if (static_rom.isNextSnapshot(0.0)) {
         static_rom.takeSnapshot(M[i]);
         double next_snapshot_time =
            static_rom.computeNextSnapshotTime(M[i], M[i], 0.0);
      }
   }
#ifdef DEBUG
   Mat static_model = static_rom.getModel();
   Mat inc_model = inc_rom.getModel();
   Mat test;
   Mat static_model_t;
   MatTranspose(static_model, MAT_INITIAL_MATRIX, &static_model_t);
   MatMatMult(static_model_t, inc_model, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &test);
   MatView(test, PETSC_VIEWER_STDOUT_WORLD);
   MatDestroy(&static_model_t);
   MatDestroy(&static_model);
   MatDestroy(&inc_model);
   MatDestroy(&test);
#endif
   for (int i = 0; i < num_snapshots; ++i) {
      delete [] M[i];
   }
   delete [] M;
   return 0;
}
