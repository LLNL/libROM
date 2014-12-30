/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: A class that reads basis vectors from a file.
 *
 *****************************************************************************/

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
         else if (time < d_time_interval_start_times[d_last_basis_idx]) {
            result = true;
         }
         else if ((d_last_basis_idx != numTimeIntervals()-1) &&
                  (d_time_interval_start_times[d_last_basis_idx+1] <= time)) {
            result = true;
         }
         return result;
      }

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
