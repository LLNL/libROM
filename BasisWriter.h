/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: A class that writes basis vectors to a file.
 *
 *****************************************************************************/

#ifndef included_BasisWriter_h
#define included_BasisWriter_h

#include "Database.h"
#include <string>

namespace CAROM {

class SVDBasisGenerator;

class BasisWriter {
   public:
      // Constructor.
      BasisWriter(
         SVDBasisGenerator* basis_generator,
         const std::string& base_file_name,
         Database::formats db_format = Database::HDF5);

      // Destructor.
      ~BasisWriter();

      // Write basis vectors for supplied SVDBasisGenerator to file with
      // specified name.
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

      // Basis generator whose basis vectors are being written.
      SVDBasisGenerator* d_basis_generator;

      // Database to which basis vectors are being written.
      Database* d_database;

      // Number of time intervals for which basis vectors have been written.
      int d_num_intervals_written;
};

}

#endif
