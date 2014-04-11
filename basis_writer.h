#ifndef basis_writer_h
#define basis_writer_h

#include "Database.h"
#include <string>

namespace CAROM {

class svd_rom;

class BasisWriter {
   public:
      // Constructor.
      BasisWriter(
         svd_rom* rom,
         const std::string& base_file_name,
         Database::formats db_format = Database::HDF5);

      // Destructor.
      ~BasisWriter();

      // Write basis vectors for supplied svd_rom to file with specified name.
      void
      writeBasis();

   private:
      // Unimplemented default constructor.
      BasisWriter();

      // Unimplemented copy constructor.
      BasisWriter(
         const BasisWriter& other);

      // Unimplemented assignment operator.
      BasisWriter&
      operator = (
         const BasisWriter& rhs);

      // ROM whose basis vectors are being written.
      svd_rom* d_rom;

      // Database to which basis vectors are being written.
      Database* d_database;

      // Number of time intervals for which basis vectors have been written.
      int d_num_intervals_written;
};

}

#endif
