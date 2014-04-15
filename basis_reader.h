#ifndef basis_reader_h
#define basis_reader_h

#include "Utilities.h"
#include "Database.h"
#include <string>
#include <vector>

namespace CAROM {

class Matrix;
class Database;

class BasisReader {
   public:
      // Constructor.
      BasisReader(
         const std::string& base_file_name,
         Database::formats db_format = Database::HDF5);

      // Destructor.
      ~BasisReader();

      // Returns true if the basis vectors at time t are different from the
      // last requested basis vectors.
      bool
      isNewBasis(
         double time)
      {
         CAROM_ASSERT(0 < d_num_time_intervals);
         CAROM_ASSERT(0 <= time);
         bool result = false;
         if (d_last_basis_idx == -1) {
            result = true;
         }
         else if ((time < d_time_interval_start_times[d_last_basis_idx])) {
            result = true;
         }
         else if ((d_last_basis_idx != d_num_time_intervals-1) &&
                  (d_time_interval_start_times[d_last_basis_idx+1] <= time)) {
            result = true;
         }
         return result;
      }

      // Basis accessor.
      const Matrix*
      getBasis(
         double time);

   private:
      // Unimplemented default constructor.
      BasisReader();

      // Unimplemented copy constructor.
      BasisReader(
         const BasisReader& other);

      // Unimplemented assignment operator.
      BasisReader&
      operator = (
         const BasisReader& rhs);

      // The number of time intervals on which basis vectors have been
      // constructed.
      int d_num_time_intervals;

      // The start time of each time interval.
      std::vector<double> d_time_interval_start_times;

      // The currently requested basis vectors.
      Matrix* d_basis_vectors;

      // The database being read from.
      Database* d_database;

      // The last time at which basis vectors were requested.
      int d_last_basis_idx;

      // The rank of the process this object belongs to.
      int d_rank;

      // The number of processors.
      int d_size;
};

}

#endif
