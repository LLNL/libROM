#ifndef svd_rom_h
#define svd_rom_h

#include "basis_writer.h"

#include <string.h>

namespace CAROM {

class BasisWriter;
class Matrix;

// An abstract base class defining the interface to the incremental svd
// algorithm.
class svd_rom
{
   public:
      // Constructor.
      svd_rom(
         const std::string& basis_file_name,
         Database::formats file_format = Database::HDF5);

      // Destructor.
      virtual
      ~svd_rom();

      // Returns true if it is time for the next svd snapshot.
      virtual
      bool
      isNextSnapshot(
         double time) = 0;

      // Add a snapshot to the incremental svd at the given time.
      virtual
      void
      takeSnapshot(
         double* u_in,
         double time) = 0;

      // Signal that the final snapshot has been taken.
      void
      endSnapshots()
      {
         if (d_basis_writer) {
            d_basis_writer->writeBasis();
         }
      }

      // Computes next time an svd snapshot is needed.
      virtual
      double
      computeNextSnapshotTime(
         double* u_in,
         double* rhs_in,
         double time) = 0;

      // Returns the basis vectors for the current time interval as a Matrix.
      virtual
      const Matrix*
      getBasis() = 0;

      // Returns the number of time intervals on which different sets of basis
      // vectors are defined.
      virtual
      int
      getNumBasisTimeIntervals() const = 0;

      // Returns the start time for the requested time interval.
      virtual
      double
      getBasisIntervalStartTime(
         int which_interval) const = 0;

   protected:
      // Writer of basis vectors.
      BasisWriter* d_basis_writer;

   private:
      // Unimplemented copy constructor.
      svd_rom(
         const svd_rom& other);

      // Unimplemented assignment operator.
      svd_rom&
      operator = (
         const svd_rom& rhs);
};

}

#endif
