#include "basis_writer.h"
#include "HDFDatabase.h"
#include "matrix.h"
#include "svd_rom.h"
#include "Utilities.h"

#include "mpi.h"

namespace CAROM {

void
BasisWriter::writeBasis(
   const std::string& base_file_name,
   svd_rom& rom,
   Database::formats db_format)
{
   CAROM_ASSERT(!base_file_name.empty());

   int rank;
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   }
   else {
      rank = 0;
   }
   char tmp[10];
   sprintf(tmp, ".%06d", rank);
   std::string full_file_name = base_file_name + tmp;
   Database* database;
   if (db_format == Database::HDF5) {
      database = new HDFDatabase();
   }
   database->create(full_file_name);
   int num_time_intervals = rom.getNumBasisTimeIntervals();
   database->putInteger("num_time_intervals", num_time_intervals);
   for (int i = 0; i < num_time_intervals; ++i) {
      double time_interval_start_time = rom.getBasisIntervalStartTime(i);
      database->putDouble("time", time_interval_start_time);
      const Matrix* basis = rom.getBasis(time_interval_start_time);
      int num_rows = basis->numRows();
      database->putInteger("num_rows", num_rows);
      int num_cols = basis->numColumns();
      database->putInteger("num_cols", num_cols);
      database->putDoubleArray(
         "basis",
         &basis->item(0, 0),
         num_rows*num_cols);
   }
   database->close();
   delete database;
}

}
