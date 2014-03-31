#include "basis_reader.h"
#include "HDFDatabase.h"
#include "matrix.h"
#include "mpi.h"

namespace CAROM {

BasisReader::BasisReader(
   const std::string& base_file_name,
   Database::formats db_format) :
   d_last_basis_idx(-1)
{
   CAROM_ASSERT(!base_file_name.empty());

   int rank;
   int size;
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
   }
   else {
      rank = 0;
      size = 1;
   }
   char tmp[10];
   sprintf(tmp, ".%06d", rank);
   std::string full_file_name = base_file_name + tmp;
   Database* database;
   if (db_format == Database::HDF5) {
      database = new HDFDatabase();
   }
   database->open(full_file_name);
   database->getInteger("num_time_intervals", d_num_time_intervals);
   d_time_interval_start_times.resize(d_num_time_intervals);
   d_basis_vectors.resize(d_num_time_intervals, 0);
   for (int i = 0; i < d_num_time_intervals; ++i) {
      database->getDouble("time", d_time_interval_start_times[i]);
      int num_rows;
      database->getInteger("num_rows", num_rows);
      int num_cols;
      database->getInteger("num_cols", num_cols);
      d_basis_vectors[i] = new Matrix(num_rows, num_cols, false, rank, size);
      database->getDoubleArray("basis",
                               &d_basis_vectors[i]->item(0, 0),
                               num_rows*num_cols);
   }
   database->close();
   delete database;
}

BasisReader::~BasisReader()
{
   for (int i = 0; i < d_num_time_intervals; ++i) {
      delete d_basis_vectors[i];
   }
}

const Matrix*
BasisReader::getBasis(
   double time)
{
   CAROM_ASSERT(0 < d_num_time_intervals);
   CAROM_ASSERT(0 <= time);
   int i;
   for (i = 0; i < d_num_time_intervals-1; ++i) {
      if (d_time_interval_start_times[i] <= time &&
          d_time_interval_start_times[i+1] < time) {
         break;
      }
   }
   d_last_basis_idx = i;
   return d_basis_vectors[i];
}

}
