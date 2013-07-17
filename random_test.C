#include "incremental_svd_rom.h"
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
//   CAROM::incremental_svd_rom rom(&argc, &argv, dim, 1.0e-12, true, 1);
   CAROM::incremental_svd_rom rom(&argc, &argv, dim, 1.0e-12, false, 1);
   int size;
   MPI_Comm_size(PETSC_COMM_WORLD, &size);
   int rank;
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   srand(1);
   for (int i = 0; i < dim*num_snapshots*rank; ++i) {
      double random = rand();
      random = random/RAND_MAX;
   }
   double** M = new double* [dim];
   for (int i = 0; i < dim; ++i) {
      M[i] = new double [num_snapshots];
      for (int j = 0; j < num_snapshots; ++j) {
         double random = rand();
         random = random/RAND_MAX;
         M[i][j] = random;
      }
   }
   for (int i = 0; i < dim*num_snapshots*(size-rank-1); ++i) {
      double random = rand();
      random = random/RAND_MAX;
   }
   for (int i = 0; i < num_lin_dep_snapshots; ++i) {
      int col = num_snapshots - i - 1;
      for (int j = 0; j < dim; ++j) {
         M[j][col] = 0;
      }
      for (int j = 0; j < num_snapshots - num_lin_dep_snapshots; ++j) {
         double random = rand();
         random = random/RAND_MAX;
         for (int k = 0; k < dim; ++k) {
            M[k][col] += random*M[k][j];
         }
      }
   }
   double* c = new double [dim];
   double* rhs = new double [dim];
   for (int i = 0; i < num_snapshots; ++i) {
      for (int j = 0; j < dim; ++j) {
         c[j] = M[j][i];
         rhs[j] = c[j];
      }
      if (rom.isNextSnapshot(0.0)) {
         rom.takeSnapshot(c);
         double next_snapshot_time =
            rom.computeNextSnapshotTime(c, rhs, 0.0);
      }
   }
   delete [] rhs;
   delete [] c;
   for (int i = 0; i < dim; ++i) {
     delete [] M[i];
   }
   delete [] M;
   return 0;
}
