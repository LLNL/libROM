/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The concrete database implementation using HDF5.
 *
 *****************************************************************************/

#ifndef included_HDFDatabase_h
#define included_HDFDatabase_h

#include "Database.h"
#include "hdf5.h"
#include <string>

namespace CAROM {

/* A class that provides basic ability to write to and read from an HDF file.
 * It's capabilities are limited to what the SVD algorithm needs to read and
 * write its basis vectors. */
class HDFDatabase : public Database
{
   public:
      // Default constructor.
      HDFDatabase();

      // Destructor.
      virtual
      ~HDFDatabase();

      // Creates a new database file with the supplied name.
      virtual
      bool
      create(
         const std::string& file_name);

      // Opens the existing database file with the supplied name.
      virtual
      bool
      open(
         const std::string& file_name);

      // Closes the database file with the supplied name.
      virtual
      bool
      close();

      // Writes an array of integers to the file.
      virtual
      void
      putIntegerArray(
         const std::string& key,
         const int* const data,
         int nelements);

      // Writes an array of doubles to the file.
      virtual
      void
      putDoubleArray(
         const std::string& key,
         const double* const data,
         int nelements);

      // Reads an array of integers from the file.
      virtual
      void
      getIntegerArray(
         const std::string& key,
         int* data,
         int nelements);

      // Reads an array of doubles from the file.
      virtual
      void
      getDoubleArray(
         const std::string& key,
         double* data,
         int nelements);

   private:
      // Unimplemented copy constructor.
      HDFDatabase(
         const HDFDatabase& other);

      // Unimplemented assignment operator.
      HDFDatabase&
      operator = (
         const HDFDatabase& rhs);

      // Returns true if the specified key represents an integer entry.  If the
      // key does not exist or if the string is empty then false is returned.
      bool
      isInteger(
         const std::string& key);

      // Returns true if the specified key represents a double entry.  If the
      // key does not exist or if the string is empty then false is returned.
      bool
      isDouble(
         const std::string& key);

      // Write attribute to the specified dataset.
      void
      writeAttribute(
         int type_key,
         hid_t dataset_id);

      // Read attribute from the specified dataset.
      int
      readAttribute(
         hid_t dataset_id);

      /* Whether database is mounted to a file. */
      bool d_is_file;

      /* ID of file attached to database or -1 if not mounted to a file. */
      hid_t d_file_id;

      /* ID of group attached to database. */
      hid_t d_group_id;

      static const int KEY_DOUBLE_ARRAY;
      static const int KEY_INT_ARRAY;
};

}

#endif
