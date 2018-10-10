/******************************************************************************
 *
 * Copyright (c) 2013-2017, Lawrence Livermore National Security, LLC.
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

#ifndef included_BasisReader_h
#define included_BasisReader_h

#include "Utilities.h"
#include "Database.h"
#include <string>
#include <vector>

namespace CAROM {

class Matrix;
class Database;

/**
 * Class BasisReader reads the basis vectors from a file written by class
 * BasisWriter.
 *
 * @see BasisWriter
 */
class BasisReader {
   public:
      /**
       * @brief The constructor for the BasisReader class takes the base part
       * of the name of the files holding the basis vectors and the file
       * format.
       *
       * @pre !base_file_name.empty()
       *
       * @param[in] base_file_name The base part of the name of the files
       *                           holding the basis vectors.
       * @param[in] db_format Format of the file to read.
       *                      One of the implemented file formats defined in
       *                      Database.
       */
      BasisReader(
         const std::string& base_file_name,
         Database::formats db_format = Database::HDF5);

      /**
       * @brief Destructor.
       */
      ~BasisReader();

      /**
       * @brief Returns true if the basis vectors at requested time are
       * different from the last requested basis vectors.
       *
       * @pre 0 < numTimeIntervals()
       * @pre 0 <= time
       *
       * @param[in] time Time at which we are interested in the basis vectors.
       *
       * @return True if the basis vectors at the requested time are different
       *         from the last requested basis vectors.
       */
      bool
      isNewBasis(
         double time)
      {
         CAROM_ASSERT(0 < numTimeIntervals());
         CAROM_ASSERT(0 <= time);
         bool result = false;
         if (d_last_basis_idx == -1) {
            result = true;
         }
         else {
            int num_time_intervals = numTimeIntervals();
            int i;
            for (i = 0; i < num_time_intervals-1; ++i) {
               if (d_time_interval_start_times[i] <= time &&
                   time < d_time_interval_start_times[i+1]) {
                  break;
               }
            }
            result = i != d_last_basis_idx;
         }
         return result;
      }

      void readBasis(
             const std::string& base_file_name,
             Database::formats db_format = Database::HDF5);

      /**
       *
       * @brief Returns the basis vectors for the requested time.
       *
       * @pre 0 < numTimeIntervals()
       * @pre 0 <= time
       *
       * @param[in] time Time for which we want the basis vectors.
       *
       * @return The basis vectors time the requested time.
       */
      const Matrix*
      getBasis(
         double time);

      const Matrix*
      getTemporalBasis(
         double time);

      Matrix
      getMatlabBasis(
         double time);

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      BasisReader();

      /**
       * @brief Unimplemented copy constructor.
       */
      BasisReader(
         const BasisReader& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      BasisReader&
      operator = (
         const BasisReader& rhs);

      /**
       * @brief Number of time intervals.
       *
       * @return The number of time intervals.
       */
      int
      numTimeIntervals()
      {
         return static_cast<int>(d_time_interval_start_times.size());
      }

      /**
       * @brief The start time of each time interval.
       */
      std::vector<double> d_time_interval_start_times;

      /**
       * @brief The currently requested basis vectors.
       */
      Matrix* d_basis_vectors;

      /**
       * @brief The currently requested temporal basis vectors.
       */
      Matrix* d_temporal_basis_vectors;

      /**
       * @brief The database being read from.
       */
      Database* d_database;

      /**
       * @brief The last time at which basis vectors were requested.
       */
      int d_last_basis_idx;
};

}

#endif
