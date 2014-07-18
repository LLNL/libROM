#ifndef included_IncrementalSVDTimeStepper_h
#define included_IncrementalSVDTimeStepper_h

#include "IncrementalSVD.h"
#include <boost/shared_ptr.hpp>

namespace CAROM {

// The class knows, given an incremental svd implementation, the time at which
// the next increment to the incremental svd is needed.  It also knows given a
// time whether it is time for the next svd increment.  There are two factors
// determining if it is time for the next svd increment:
// 1) the current time compared to the time the next increment must happen
// 2) the number of time steps since the last increment
class IncrementalSVDTimeStepper
{
   public:
      // Constructor.
      IncrementalSVDTimeStepper(
         int dim,
         double redundancy_tol,
         bool skip_redundant,
         int increments_per_time_interval,
         double sampling_tol,
         double max_time_between_increments,
         bool fast_update);

      // Destructor.
      ~IncrementalSVDTimeStepper();

      // Returns true if it is time for the next svd increment.
      bool
      isNextIncrement(
         double time)
      {
         return time >= d_next_increment_time;
      }

      // Increment the incremental svd at the given time.
      void
      increment(
         double* u_in,
         double time)
      {
         d_isvd->increment(u_in, time);
      }

      // Computes next time an svd increment is needed.
      double
      computeNextIncrementTime(
         double* u_in,
         double* rhs_in,
         double time);

      // Returns the basis vectors for the current time interval as a Matrix.
      const Matrix*
      getBasis()
      {
         return d_isvd->getBasis();
      }

      // Returns the number of time intervals on which different sets of basis
      // vectors are defined.
      int
      getNumBasisTimeIntervals() const
      {
         return d_isvd->getNumBasisTimeIntervals();
      }

      // Returns the start time for the requested time interval.
      double
      getBasisIntervalStartTime(
         int which_interval) const
      {
         return d_isvd->getBasisIntervalStartTime(which_interval);
      }

      // Returns true if the next state will result in a new time interval.
      bool
      isNewTimeInterval() const
      {
         return d_isvd->isNewTimeInterval();
      }

   private:
      // Unimplemented default constructor.
      IncrementalSVDTimeStepper();

      // Unimplemented copy constructor.
      IncrementalSVDTimeStepper(
         const IncrementalSVDTimeStepper& other);

      // Unimplemented assignment operator.
      IncrementalSVDTimeStepper&
      operator = (
         const IncrementalSVDTimeStepper& rhs);

      // Tolerance for norm of error in the solution of the reduced model due
      // to orthogonal projection.
      double d_tol;

      // Maximum time between increments.
      double d_max_time_between_increments;

      // Next time at which an increment should be taken.
      double d_next_increment_time;

      // The fundamental incremental SVD algorithm.
      boost::shared_ptr<IncrementalSVD> d_isvd;
};

}

#endif
