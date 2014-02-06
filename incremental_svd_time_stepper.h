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
         double tolerance,
         double max_time_between_increments,
         bool fast_update);

      // Destructor.
      ~incremental_svd_time_stepper();

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

      // Returns the model parameters at the given time as a Matrix.
      const Matrix*
      getModel(
         double time)
      {
         return d_isvd->getModel(time);
      }

      // Writes the model to a file with the given base name.
      void
      writeModel(
         const std::string& base_file_name)
      {
         d_isvd->writeModel(base_file_name);
      }

      // Reads the model from a file with the given base name.
      void
      readModel(
         const std::string& base_file_name)
      {
         d_isvd->readModel(base_file_name);
      }

      double
      checkOrthogonality()
      {
         return d_isvd->checkOrthogonality();
      }

      void
      reOrthogonalize()
      {
         d_isvd->reOrthogonalize();
      }

   private:
      // Unimplemented default constructor.
      incremental_svd_time_stepper();

      // Unimplemented copy constructor.
      incremental_svd_time_stepper(
         const incremental_svd_time_stepper& other);

      // Unimplemented assignment operator.
      incremental_svd_time_stepper&
      operator = (
         const incremental_svd_time_stepper& rhs);

      // Tolerance for norm of error in the solution of the reduced model due
      // to orthogonal projection.
      double d_tol;

      // Maximum time between increments.
      double d_max_time_between_increments;

      // Next time at which an increment should be taken.
      double d_next_increment_time;

      // The fundamental incremental SVD algorithm.
      boost::shared_ptr<incremental_svd> d_isvd;
};

}

#endif
