/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2015 Lawrence Livermore National Security, LLC
 * Description: A class that writes basis vectors to a file.
 *
 *****************************************************************************/

#include "BasisWriter.h"
#include "HDFDatabase.h"
#include "Matrix.h"
#include "SVDBasisGenerator.h"
#include "Utilities.h"

#include "mpi.h"

namespace CAROM {

BasisWriter::BasisWriter(
   SVDBasisGenerator* basis_generator,
   const std::string& base_file_name,
   Database::formats db_format) :
   d_basis_generator(basis_generator),
   d_num_intervals_written(0)
{
   CAROM_ASSERT(basis_generator != 0);
   CAROM_ASSERT(!base_file_name.empty());

   int mpi_init;
   MPI_Initialized(&mpi_init);
   int rank;
   if (mpi_init) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   }
   else {
      rank = 0;
   }
   char tmp[100];
   sprintf(tmp, ".%06d", rank);
   std::string full_file_name = base_file_name + tmp;
   if (db_format == Database::HDF5) {
      d_database = new HDFDatabase();
   }
   d_database->create(full_file_name);
}

BasisWriter::~BasisWriter()
{
   d_database->putInteger("num_time_intervals", d_num_intervals_written);
   d_database->close();
   delete d_database;
}

void
BasisWriter::writeBasis()
{
   char tmp[100];
   double time_interval_start_time =
      d_basis_generator->getBasisIntervalStartTime(d_num_intervals_written);
   sprintf(tmp, "time_%06d", d_num_intervals_written);
   d_database->putDouble(tmp, time_interval_start_time);
   const Matrix* basis = d_basis_generator->getBasis();
   int num_rows = basis->numRows();
   sprintf(tmp, "num_rows_%06d", d_num_intervals_written);
   d_database->putInteger(tmp, num_rows);
   int num_cols = basis->numColumns();
   sprintf(tmp, "num_cols_%06d", d_num_intervals_written);
   d_database->putInteger(tmp, num_cols);
   sprintf(tmp, "basis_%06d", d_num_intervals_written);
   d_database->putDoubleArray(tmp, &basis->item(0, 0), num_rows*num_cols);
   ++d_num_intervals_written;
}

}
