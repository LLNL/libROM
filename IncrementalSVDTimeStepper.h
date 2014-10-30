/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The class that determines the next time at which a sample
 *              should be taken for basis generation using an incremental SVD
 *              approach.
 *
 *****************************************************************************/

#ifndef included_IncrementalSVDTimeStepper_h
#define included_IncrementalSVDTimeStepper_h

#include "IncrementalSVD.h"
#include <boost/shared_ptr.hpp>

namespace CAROM {

/**
 * Class IncrementalSVDTimeStepper knows, given an incremental svd
 * implementation, the time at which the next increment to the incremental svd
 * is needed.  It also knows given a time whether it is time for the next svd
 * increment.  There are two factors determining if it is time for the next
 * svd increment:
 * 1) the current time compared to the time the next increment must happen
 * 2) the number of time steps since the last increment
 */
class IncrementalSVDTimeStepper
{
   public:
      /**
       * @brief Constructor.
       *
       * @pre dim > 0
       * @pre redundancy_tol > 0.0
       * @pre increments_per_time_interval > 0
       * @pre sampling_tol > 0.0
       * @pre max_time_between_increments > 0.0
       *
       * @param[in] dim The dimension of the system on this processor.
       * @param[in] redundancy_tol Tolerance to determine if a sample is
       *                           redundant or not.
       * @param[in] skip_redundant If true skip redundant samples.
       * @param[in] increments_per_time_interval The maximum number of samples
       *                                         in each time interval.
       * @param[in] sampling_tol Time step control tolerance.  Limits error in
       *                         projection of solution into reduced order
       *                         space.
       * @param[in] max_time_between_increments Hard upper bound on time step.
       * @param[in] fast_update If true use the fast update incremental svd
       *                        algorithm.
       */
      IncrementalSVDTimeStepper(
         int dim,
         double redundancy_tol,
         bool skip_redundant,
         int increments_per_time_interval,
         double sampling_tol,
         double max_time_between_increments,
         bool fast_update);

      /**
       * @brief Destructor.
       */
      ~IncrementalSVDTimeStepper();

      /**
       * @brief Returns true if it is time for the next svd increment.
       *
       * @param[in] time Time of interest.
       *
       * @return True if it is time for the next snapshot to be taken.
       */
      bool
      isNextIncrement(
         double time)
      {
         return time >= d_next_increment_time;
      }

      /**
       * @brief Increment the incremental svd at the given time with the
       * supplied state.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] time The simulation time for the state.
       */
      void
      increment(
         double* u_in,
         double time)
      {
         d_isvd->increment(u_in, time);
      }

      /**
       * @brief Computes next time an svd increment is needed.
       *
       * @pre u_in != 0
       * @pre rhs_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] rhs_in The right hand side at the specified time.
       * @param[in] time The simulation time for the state.
       */
      double
      computeNextIncrementTime(
         double* u_in,
         double* rhs_in,
         double time);

      /**
       * @brief Returns the basis vectors for the current time interval as a
       * Matrix.
       *
       * @return The basis vectors for the current time interval.
       */
      const Matrix*
      getBasis()
      {
         return d_isvd->getBasis();
      }

      /**
       * @brief Returns the number of time intervals on which different sets of
       * basis vectors are defined.       *
       * @return The number of time intervals on which there are basis vectors.
       */
      int
      getNumBasisTimeIntervals() const
      {
         return d_isvd->getNumBasisTimeIntervals();
      }

      /**
       * @brief Returns the start time for the requested time interval.
       *
       * @pre 0 <= which_interval
       * @pre which_interval < getNumBasisTimeIntervals()
       *
       * @param[in] which_interval Time interval whose start time is needed.
       *
       * @return The start time for the requested time interval.
       */
      double
      getBasisIntervalStartTime(
         int which_interval) const
      {
         return d_isvd->getBasisIntervalStartTime(which_interval);
      }

      /**
       * @brief Returns true if the next state will result in a new time
       * interval.
       *
       * @return True if all samples have been taken for the current time
       * interval.
       */
      bool
      isNewTimeInterval() const
      {
         return d_isvd->isNewTimeInterval();
      }

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      IncrementalSVDTimeStepper();

      /**
       * @brief Unimplemented copy constructor.
       */
      IncrementalSVDTimeStepper(
         const IncrementalSVDTimeStepper& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      IncrementalSVDTimeStepper&
      operator = (
         const IncrementalSVDTimeStepper& rhs);

      /**
       * @brief Time step control tolerance.
       *
       * Limits error in projection of solution into the reduced order space.
       */
      double d_tol;

      /**
       * @brief Maximum time between increments.
       */
      double d_max_time_between_increments;

      /**
       * @brief Next time at which an increment should be taken.
       */
      double d_next_increment_time;

      /**
       * @brief Pointer to the fundamental incremental SVD algorithm object.
       */
      boost::shared_ptr<IncrementalSVD> d_isvd;
};

}

#endif
