/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The abstract wrapper class for an abstract SVD algorithm and
 *              time stepper.  This class provides interfaces to each so that
 *              an application only needs to instantiate one concrete
 *              implementation of this class to control all aspects of basis
 *              vector generation.
 *
 *****************************************************************************/

#ifndef included_SVDBasisGenerator_h
#define included_SVDBasisGenerator_h

#include "BasisWriter.h"

#include <string.h>

namespace CAROM {

class BasisWriter;
class Matrix;

// An abstract base class defining the interface to the incremental svd
// algorithm.
class SVDBasisGenerator
{
   public:
      // Constructor.
      SVDBasisGenerator(
         const std::string& basis_file_name,
         Database::formats file_format = Database::HDF5);

      // Destructor.
      virtual
      ~SVDBasisGenerator();

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
      SVDBasisGenerator(
         const SVDBasisGenerator& other);

      // Unimplemented assignment operator.
      SVDBasisGenerator&
      operator = (
         const SVDBasisGenerator& rhs);
};

}

#endif
