#ifndef included_static_svd_basis_generator_h
#define included_static_svd_basis_generator_h

#include "svd_basis_generator.h"
#include "static_svd_time_stepper.h"

namespace CAROM {

// A class which implements a Reduced Order Model based on the static svd
// algorithm.
class static_svd_basis_generator : public svd_basis_generator
{
   public:
      // Constructor.
      static_svd_basis_generator(
         int dim,
         int increments_per_time_interval,
         const std::string& basis_file_name,
         Database::formats file_format = Database::HDF5);

      // Destructor.
      virtual
      ~static_svd_basis_generator();

      // Returns true if it is time for the next svd snapshot.
      virtual
      bool
      isNextSnapshot(
         double time);

      // Add a snapshot to the static svd.
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
      static_svd_basis_generator();

      // Unimplemented copy constructor.
      static_svd_basis_generator(
         const static_svd_basis_generator& other);

      // Unimplemented assignment operator.
      static_svd_basis_generator&
      operator = (
         const static_svd_basis_generator& rhs);

      boost::shared_ptr<static_svd_time_stepper> d_svdts;
};

}

#endif
