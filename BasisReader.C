/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: A class that reads basis vectors from a file.
 *
 *****************************************************************************/

#include "BasisReader.h"
#include "HDFDatabase.h"
#include "Matrix.h"
#include "mpi.h"

namespace CAROM {

BasisReader::BasisReader(
   const std::string& base_file_name,
   Database::formats db_format) :
   d_basis_vectors(0),
   d_last_basis_idx(-1)
{
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
   d_database->open(full_file_name);
   int num_time_intervals;
   d_database->getInteger("num_time_intervals", num_time_intervals);
   d_time_interval_start_times.resize(num_time_intervals);
   for (int i = 0; i < num_time_intervals; ++i) {
      sprintf(tmp, "time_%06d", i);
      d_database->getDouble(tmp, d_time_interval_start_times[i]);
   }
}

BasisReader::~BasisReader()
{
   if (d_basis_vectors) {
      delete d_basis_vectors;
   }
   d_database->close();
   delete d_database;
}

const Matrix*
BasisReader::getBasis(
   double time)
{
   CAROM_ASSERT(0 < numTimeIntervals());
   CAROM_ASSERT(0 <= time);
   int num_time_intervals = numTimeIntervals();
   int i;
   for (i = 0; i < num_time_intervals-1; ++i) {
      if (d_time_interval_start_times[i] <= time &&
          time < d_time_interval_start_times[i+1]) {
         break;
      }
   }
   d_last_basis_idx = i;
   char tmp[100];
   int num_rows;
   sprintf(tmp, "num_rows_%06d", i);
   d_database->getInteger(tmp, num_rows);
   int num_cols;
   sprintf(tmp, "num_cols_%06d", i);
   d_database->getInteger(tmp, num_cols);
   if (d_basis_vectors) {
      delete d_basis_vectors;
   }
   d_basis_vectors = new Matrix(num_rows, num_cols, true);
   sprintf(tmp, "basis_%06d", i);
   d_database->getDoubleArray(tmp,
                              &d_basis_vectors->item(0, 0),
                              num_rows*num_cols);
   return d_basis_vectors;
}

}
