/******************************************************************************
 *
 * Copyright (c) 2013-2016, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 * Written by William Arrighi wjarrighi@llnl.gov
 * CODE-686965
 * All rights reserved.
 *
 * This file is part of libROM.
 * For details, see https://computation.llnl.gov/librom
 * Please also read README_BSD_NOTICE.
 *
 * Redistribution and use in source and binary forms, with or without
 * modifications, are permitted provided that the following conditions are met:
 *
 *    o Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the disclaimer below.
 *    o Redistribution in binary form must reproduce the above copyright
 *      notice, this list of conditions and the disclaimer (as noted below) in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *    o Neither the name of the LLNS/LLNL nor the names of its contributors may
 *      be used to endorse or promote products derived from this software
 *      without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
 * LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OR SUCH DAMAGE.
 *
 *****************************************************************************/

// Description: A class that reads basis vectors from a file.

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
