#ifndef included_IncrementalSVDBasisGenerator_h
#define included_IncrementalSVDBasisGenerator_h

#include "SVDBasisGenerator.h"
#include "IncrementalSVDTimeStepper.h"

namespace CAROM {

// A class which implements a Reduced Order Model based on the incremental svd
// algorithm.
  class IncrementalSVDBasisGenerator : public SVDBasisGenerator
{
   public:
      // Constructor.
      IncrementalSVDBasisGenerator(
         int dim,
         double redundancy_tol,
         bool skip_redundant,
         int basis_size,
         double sampling_tol,
         double max_time_between_snapshots,
         bool fast_update,
         const std::string& basis_file_name,
         Database::formats file_format = Database::HDF5);

      // Destructor.
      virtual
      ~IncrementalSVDBasisGenerator();

      // Returns true if it is time for the next svd snapshot.
      virtual
      bool
      isNextSnapshot(
         double time);

      // Add a snapshot to the incremental svd at the given time.
      virtual
      void
      takeSnapshot(
         double* u_in,
         double time);

      // Computes next time an svd snapshot is needed.
      virtual
      double
      computeNextSnapshotTime(
         double* u_in,
         double* rhs_in,
         double time);

      // Returns the basis vectors for the current time interval as a Matrix.
      virtual
      const Matrix*
      getBasis();

      // Returns the number of time intervals on which different sets of basis
      // vectors are defined.
      virtual
      int
      getNumBasisTimeIntervals() const;

      // Returns the start time for the requested time interval.
      virtual
      double
      getBasisIntervalStartTime(
         int which_interval) const;

   private:
      // Unimplemented default constructor.
      IncrementalSVDBasisGenerator();

      // Unimplemented copy constructor.
      IncrementalSVDBasisGenerator(
         const IncrementalSVDBasisGenerator& other);

      // Unimplemented assignment operator.
      IncrementalSVDBasisGenerator&
      operator = (
         const IncrementalSVDBasisGenerator& rhs);

      // The underlying time step control object.
      boost::shared_ptr<IncrementalSVDTimeStepper> d_isvdts;
};

}

#endif
