#ifndef included_Database_h
#define included_Database_h

#include <string>

namespace CAROM {

/* An abstract base class that provides basic ability to write to and read
 * from a file.  It's capabilities are limited to what the SVD algorithm needs
 * to read and write its basis vectors. */
class Database
{
   public:
      // Default constructor.
      Database();

      // Destructor.
      virtual
      ~Database();

      // Creates a new database file with the supplied name.
      virtual
      bool
      create(
         const std::string& file_name) = 0;

      // Opens the existing database file with the supplied name.
      virtual
      bool
      open(
         const std::string& file_name) = 0;

      // Closes the database file with the supplied name.
      virtual
      bool
      close() = 0;

      // Writes an integer to the file.
      void
      putInteger(
         const std::string& key,
         int data)
      {
         putIntegerArray(key, &data, 1);
      }

      // Writes an array of integers to the file.
      virtual
      void
      putIntegerArray(
         const std::string& key,
         const int* const data,
         int nelements) = 0;

      // Writes a double to the file.
      void
      putDouble(
         const std::string& key,
         double data)
      {
         putDoubleArray(key, &data, 1);
      }

      // Writes an array of doubles to the file.
      virtual
      void
      putDoubleArray(
         const std::string& key,
         const double* const data,
         int nelements) = 0;

      // Reads an integer from the file.
      void
      getInteger(
         const std::string& key,
         int& data)
      {
         getIntegerArray(key, &data, 1);
      }

      // Reads an array of integers from the file.
      virtual
      void
      getIntegerArray(
         const std::string& key,
         int* data,
         int nelements) = 0;

      // Reads a double from the file.
      void
      getDouble(
         const std::string& key,
         double& data)
      {
         getDoubleArray(key, &data, 1);
      }

      // Reads an array of doubles from the file.
      virtual
      void
      getDoubleArray(
         const std::string& key,
         double* data,
         int nelements) = 0;

      enum formats {
         HDF5
      };

   private:
      // Unimplemented copy constructor.
      Database(
         const Database& other);

      // Unimplemented assignment operator.
      Database&
      operator = (
         const Database& rhs);
};

}

#endif
