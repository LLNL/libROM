#ifndef included_incremental_svd_basis_generator_h
#define included_incremental_svd_basis_generator_h

#include "svd_basis_generator.h"
#include "incremental_svd_time_stepper.h"

namespace CAROM {

// A class which implements a Reduced Order Model based on the incremental svd
// algorithm.
  class incremental_svd_basis_generator : public svd_basis_generator
{
   public:
      // Constructor.
      incremental_svd_basis_generator(
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
      ~incremental_svd_basis_generator();

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
      incremental_svd_basis_generator();

      // Unimplemented copy constructor.
      incremental_svd_basis_generator(
         const incremental_svd_basis_generator& other);

      // Unimplemented assignment operator.
      incremental_svd_basis_generator&
      operator = (
         const incremental_svd_basis_generator& rhs);

      // The underlying time step control object.
      boost::shared_ptr<incremental_svd_time_stepper> d_isvdts;
};

}

#endif
