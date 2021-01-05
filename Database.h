/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The abstract database class defines interface for databases.

#ifndef included_Database_h
#define included_Database_h

#include <string>

namespace CAROM {

/**
 * Class Database is an abstract base class that provides basic ability to
 * write to and read from a file.  It's capabilities are limited to what the
 * SVD algorithm needs to read and write its basis vectors.
 */
class Database
{
   public:
      /**
       * @brief Default constructor.
       */
      Database();

      /**
       * @brief Destructor.
       */
      virtual
      ~Database();

      /**
       * @brief Creates a new database file with the supplied name.
       *
       * @param[in] file_name Name of database file to create.
       *
       * @return true if file create was successful.
       */
      virtual
      bool
      create(
         const std::string& file_name) = 0;

      /**
       * @brief Opens an existing database file with the supplied name.
       *
       * @param[in] file_name Name of existing database file to open.
       *
       * @return true if file open was successful.
       */
      virtual
      bool
      open(
         const std::string& file_name) = 0;

      /**
       * @brief Closes the currently open database file.
       *
       * @return true if the file close was successful.
       */
      virtual
      bool
      close() = 0;

      /**
       * @brief Writes an integer associated with the supplied key to currently
       * open database file.
       *
       * @param[in] key The key associated with the value to be written.
       * @param[in] data The integer value to be written.
       */
      void
      putInteger(
         const std::string& key,
         int data)
      {
         putIntegerArray(key, &data, 1);
      }

      /**
       * @brief Writes an array of integers associated with the supplied key to
       * the currently open database file.
       *
       * @param[in] key The key associated with the array of values to be
       *                written.
       * @param[in] data The array of integer values to be written.
       * @param[in] nelements The number of integers in the array.
       */
      virtual
      void
      putIntegerArray(
         const std::string& key,
         const int* const data,
         int nelements) = 0;

      /**
       * @brief Writes a double associated with the supplied key to currently
       * open database file.
       *
       * @param[in] key The key associated with the value to be written.
       * @param[in] data The double value to be written.
       */
      void
      putDouble(
         const std::string& key,
         double data)
      {
         putDoubleArray(key, &data, 1);
      }

      /**
       * @brief Writes an array of doubles associated with the supplied key to
       * the currently open database file.
       *
       * @param[in] key The key associated with the array of values to be
       *                written.
       * @param[in] data The array of double values to be written.
       * @param[in] nelements The number of doubles in the array.
       */
      virtual
      void
      putDoubleArray(
         const std::string& key,
         const double* const data,
         int nelements) = 0;

      /**
       * @brief Reads an integer associated with the supplied key from the
       * currently open database file.
       *
       * @param[in] key The key associated with the value to be read.
       * @param[out] data The integer value read.
       */
      void
      getInteger(
         const std::string& key,
         int& data)
      {
         getIntegerArray(key, &data, 1);
      }

      /**
       * @brief Reads an array of integers associated with the supplied key
       * from the currently open database file.
       *
       * @param[in] key The key associated with the array of values to be
       *                read.
       * @param[out] data The allocated array of integer values to be read.
       * @param[in] nelements The number of integers in the array.
       */
      virtual
      void
      getIntegerArray(
         const std::string& key,
         int* data,
         int nelements) = 0;

      /**
       * @brief Reads a double associated with the supplied key from the
       * currently open database file.
       *
       * @param[in] key The key associated with the value to be read.
       * @param[out] data The double value read.
       */
      void
      getDouble(
         const std::string& key,
         double& data)
      {
         getDoubleArray(key, &data, 1);
      }

      /**
       * @brief Reads an array of doubles associated with the supplied key
       * from the currently open database file.
       *
       * @param[in] key The key associated with the array of values to be
       *                read.
       * @param[out] data The allocated array of double values to be read.
       * @param[in] nelements The number of doubles in the array.
       */
      virtual
      void
      getDoubleArray(
         const std::string& key,
         double* data,
         int nelements) = 0;

      /**
       * @brief Implemented database file formats.
       *
       * Add to this enum each time a new database format is implemented.
       */
      enum formats {
         HDF5
      };

   private:
      /**
       * @brief Unimplemented copy constructor.
       */
      Database(
         const Database& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      Database&
      operator = (
         const Database& rhs);
};

}

#endif
