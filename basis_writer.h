#ifndef basis_writer_h
#define basis_writer_h

#include "Database.h"
#include <string>

namespace CAROM {

class svd_rom;

class BasisWriter {
   public:
      // Destructor.
      ~BasisWriter();

      // Write basis vectors for supplied svd_rom to file with specified name.
      static
      void
      writeBasis(
         const std::string& base_file_name,
         svd_rom& rom,
         Database::formats db_format = Database::HDF5);

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
};

}

#endif
