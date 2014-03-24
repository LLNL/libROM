#ifndef basis_reader_h
#define basis_reader_h

#include "Database.h"
#include <string>
#include <vector>

namespace CAROM {

class Matrix;

class BasisReader {
   public:
      // Constructor.
      BasisReader(
         const std::string& base_file_name,
         Database::formats db_format = Database::HDF5);

      // Destructor.
      ~BasisReader();

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

      // The basis vectors for each time interval.
      std::vector<Matrix*> d_basis_vectors;
};

}

#endif
