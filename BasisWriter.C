/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: A class that writes basis vectors to a file.

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

   // create and open basis database
   char tmp[100];
   sprintf(tmp, ".%06d", rank);
   std::string full_file_name = base_file_name + tmp;
   if (db_format == Database::HDF5) {
      d_database = new HDFDatabase();
   }
   d_database->create(full_file_name);
   
   // create and open state database
   char tmp2[100];
   sprintf(tmp2, ".state.%06d", rank);
   std::string state_file_name = base_file_name + tmp2;
   if (db_format == Database::HDF5) {
      d_state_database = new HDF5Database();
   }
   d_state_database->create(state_file_name);
}

BasisWriter::~BasisWriter()
{
   d_database->putInteger("num_time_intervals", d_num_intervals_written);
   d_database->close();
   d_state_database->close();
   delete d_database;
   delete d_state_database;
}

void
BasisWriter::writeBasis(const std::string& kind)
{

   char tmp[100];
   double time_interval_start_time =
      d_basis_generator->getBasisIntervalStartTime(d_num_intervals_written);
   sprintf(tmp, "time_%06d", d_num_intervals_written);
   
   if (kind == "basis") {
   d_database->putDouble(tmp, time_interval_start_time);

   const Matrix* basis = d_basis_generator->getSpatialBasis();
   int num_rows = basis->numRows();
   sprintf(tmp, "spatial_basis_num_rows_%06d", d_num_intervals_written);
   d_database->putInteger(tmp, num_rows);
   int num_cols = basis->numColumns();
   sprintf(tmp, "spatial_basis_num_cols_%06d", d_num_intervals_written);
   d_database->putInteger(tmp, num_cols);
   sprintf(tmp, "spatial_basis_%06d", d_num_intervals_written);
   d_database->putDoubleArray(tmp, &basis->item(0, 0), num_rows*num_cols);

   if(d_basis_generator->updateRightSV()) {
     const Matrix* tbasis = d_basis_generator->getTemporalBasis();
     num_rows = tbasis->numRows();
     sprintf(tmp, "temporal_basis_num_rows_%06d", d_num_intervals_written);
     d_database->putInteger(tmp, num_rows);
     num_cols = tbasis->numColumns();
     sprintf(tmp, "temporal_basis_num_cols_%06d", d_num_intervals_written);
     d_database->putInteger(tmp, num_cols);
     sprintf(tmp, "temporal_basis_%06d", d_num_intervals_written);
     d_database->putDoubleArray(tmp, &tbasis->item(0, 0), num_rows*num_cols);
   }

   const Matrix* sv = d_basis_generator->getSingularValues();
   num_rows = sv->numRows();
   num_cols = num_rows;
   sprintf(tmp, "singular_value_size_%06d", d_num_intervals_written);
   d_database->putInteger(tmp, num_rows);
   sprintf(tmp, "singular_value_%06d", d_num_intervals_written);
   d_database->putDoubleArray(tmp, &sv->item(0, 0), num_rows*num_cols);

   ++d_num_intervals_written;
   }
   if (kind == "state") {
   // TODO: Implement state writing (see incremental svd code)
      

   }
   if (kind != "basis" || kind != "state") { 
      std::cout << "Basis Writer needs to write a basis or state, please call BasisWriter::writeBasis correctly"
      << std::endl;
   }
 

}

}
