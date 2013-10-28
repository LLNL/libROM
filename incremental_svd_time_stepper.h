#ifndef included_incremental_svd_time_stepper_h
#define included_incremental_svd_time_stepper_h

#include "incremental_svd.h"
#include <boost/shared_ptr.hpp>

namespace CAROM {

// The class knows, given an incremental svd implementation, the time at which
// the next increment to the incremental svd is needed.  It also knows given a
// time whether it is time for the next svd increment.  There are two factors
// determining if it is time for the next svd increment:
// 1) the current time compared to the time the next increment must happen
// 2) the number of time steps since the last increment
class incremental_svd_time_stepper
{
   public:
      // Constructor.
      incremental_svd_time_stepper(
         int dim,
         double epsilon,
         bool skip_redundant,
         int increments_per_time_interval,
         int max_time_steps_between_increments);

      // Destructor.
      ~incremental_svd_time_stepper();

      // Returns true if it is time for the next svd increment.
      bool
      isNextIncrement(
         double time)
      {
         bool result = (time >= d_next_increment_time) ||
            (d_time_steps_since_last_increment >= d_max_time_steps_between_increments);
         if (result) {
            d_time_steps_since_last_increment = 1;
         }
         else {
            ++d_time_steps_since_last_increment;
         }
         return result;
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

      // Returns the model parameters at the given time as a Matrix.
      Matrix*
      getModel(
         double time)
      {
         return d_isvd->getModel(time);
      }

   private:
      // Unimplemented default constructor.
      incremental_svd_time_stepper();

      // Unimplemented copy constructor.
      incremental_svd_time_stepper(const incremental_svd_time_stepper& other);

      // Unimplemented assignment operator.
      incremental_svd_time_stepper&
      operator = (const incremental_svd_time_stepper& rhs);

      // Maximum number time steps between increments.
      int d_max_time_steps_between_increments;

      // Number of time steps since last increment.
      int d_time_steps_since_last_increment;

      // Next time at which an increment should be taken.
      double d_next_increment_time;

      // The fundamental incremental SVD algorithm.
      boost::shared_ptr<incremental_svd> d_isvd;
};

}

#endif
