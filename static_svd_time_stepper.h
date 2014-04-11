#ifndef included_static_svd_time_stepper_h
#define included_static_svd_time_stepper_h

#include "static_svd.h"
#include <boost/shared_ptr.hpp>

namespace CAROM {

// The class knows, given a static svd implementation, the time at which
// the next state collection is needed.  It also knows given a time whether
// it is time for the next state collection.
class static_svd_time_stepper
{
   public:
      // Constructor.
      static_svd_time_stepper(
         int dim,
         int increments_per_time_interval);

      // Destructor.
      ~static_svd_time_stepper();

      // Returns true if it is time for the next state collection.
      bool
      isNextStateCollection(
         double time)
      {
         NULL_USE(time);
         return true;
      }

      // Collect the new state, u_in.
      void
      collectState(
         double* u_in,
         double time)
      {
         d_svd->collectState(u_in, time);
      }

      // Computes next time a state collection is needed.
      double
      computeNextStateCollectionTime(
         double* u_in,
         double* rhs_in,
         double time)
      {
         NULL_USE(u_in);
         NULL_USE(rhs_in);
         return time;
      }

      // Returns the basis vectors for the current time interval as a Matrix.
      const Matrix*
      getBasis()
      {
         return d_svd->getBasis();
      }

      // Returns the number of time intervals on which different sets of basis
      // vectors are defined.
      int
      getNumBasisTimeIntervals() const
      {
         return d_svd->getNumBasisTimeIntervals();
      }

      // Returns the start time for the requested time interval.
      double
      getBasisIntervalStartTime(
         int which_interval) const
      {
         return d_svd->getBasisIntervalStartTime(which_interval);
      }

      // Returns true if the next state will result in a new time interval.
      bool
      isNewTimeInterval() const
      {
         return d_svd->isNewTimeInterval();
      }

   private:
      // Unimplemented default constructor.
      static_svd_time_stepper();

      // Unimplemented copy constructor.
      static_svd_time_stepper(
         const static_svd_time_stepper& other);

      // Unimplemented assignment operator.
      static_svd_time_stepper&
      operator = (
         const static_svd_time_stepper& rhs);

      // The fundamental static SVD algorithm.
      boost::shared_ptr<static_svd> d_svd;
};

}

#endif
