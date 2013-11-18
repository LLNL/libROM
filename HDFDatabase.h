#ifndef included_HDFDatabase_h
#define included_HDFDatabase_h

#include "hdf5.h"
#include <string>

namespace CAROM {

/* A class that provides basic ability to write to and read from an HDF file.
 * It's capabilities are limited to what the incremental SVD algorithm needs
 * to read and write its model data. */
class HDFDatabase
{
   public:
      // Default constructor.
      HDFDatabase();

      // Destructor.
      ~HDFDatabase();

      bool
      create(
         const std::string& file_name);

      bool
      open(
         const std::string& file_name);

      bool
      close();

      // Writes an integer to the file.
      void
      putInteger(
         const std::string& key,
         int data)
      {
         putIntegerArray(key, &data, 1);
      }

      // Writes an array of integers to the file.
      void
      putIntegerArray(
         const std::string& key,
         const int* const data,
         int nelements);

      // Writes a double to the file.
      void
      putDouble(
         const std::string& key,
         double data)
      {
         putDoubleArray(key, &data, 1);
      }

      // Writes an array of doubles to the file.
      void
      putDoubleArray(
         const std::string& key,
         const double* const data,
         int nelements);

      // Reads an integer from the file.
      void
      getInteger(
         const std::string& key,
         int& data)
      {
         getIntegerArray(key, &data, 1);
      }

      // Reads an array of integers from the file.
      void
      getIntegerArray(
         const std::string& key,
         int* data,
         int nelements);

      // Reads a double from the file.
      void
      getDouble(
         const std::string& key,
         double& data)
      {
         getDoubleArray(key, &data, 1);
      }

      // Reads an array of doubles from the file.
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

      bool
      isInteger(
         const std::string& key);

      bool
      isDouble(
         const std::string& key);

      void
      writeAttribute(
         int type_key,
         hid_t dataset_id);

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
