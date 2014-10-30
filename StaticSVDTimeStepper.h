/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The class that determines the next time at which a sample
 *              should be taken for basis generation using the static SVD
 *              approach.
 *
 *****************************************************************************/

#ifndef included_StaticSVDTimeStepper_h
#define included_StaticSVDTimeStepper_h

#include "StaticSVD.h"
#include <boost/shared_ptr.hpp>

namespace CAROM {

/**
 * Class StaticSVDTimeStepper knows, given a static svd implementation, the
 * time at which the next state collection is needed.  It also knows given a
 * time whether it is time for the next state collection.  All states are
 * sampled in the static SVD implementation so it is always time for a new
 * state collection.
 */
class StaticSVDTimeStepper
{
   public:
      /**
       * @brief Constructor.
       *
       * @pre dim > 0
       * @pre increments_per_time_interval > 0
       *
       * @param[in] dim The dimension of the system on this processor.
       * @param[in] increments_per_time_interval The maximum number of samples
       *                                         in each time interval.
       */
      StaticSVDTimeStepper(
         int dim,
         int increments_per_time_interval);

      /**
       * @brief Destructor.
       */
      ~StaticSVDTimeStepper();

      /**
       * @brief Returns true if it is time for the next state collection.
       *
       * As the static algorithm samples everything this always returns true.
       *
       * @param[in] time Time of interest--unused.
       *
       * @return true
       */
      bool
      isNextStateCollection(
         double time)
      {
         CAROM_NULL_USE(time);
         return true;
      }

      /**
       * @brief Collect the new state, u_in.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] time The simulation time for the state.
       */
      void
      collectState(
         double* u_in,
         double time)
      {
         d_svd->collectState(u_in, time);
      }

      /**
       * @brief Computes next time a state collection is needed.
       *
       * @param[in] u_in The state at the specified time--unused.
       * @param[in] rhs_in The right hand side at the specified time--unused.
       * @param[in] time The simulation time for the state.
       *
       * @return The current simulation time as the static algorithm samples at
       * each time step.
       */
      double
      computeNextStateCollectionTime(
         double* u_in,
         double* rhs_in,
         double time)
      {
         CAROM_NULL_USE(u_in);
         CAROM_NULL_USE(rhs_in);
         return time;
      }

      /**
       * @brief Returns the basis vectors for the current time interval as a
       * Matrix.
       *
       * @return The basis vectors for the current time interval.
       */
      const Matrix*
      getBasis()
      {
         return d_svd->getBasis();
      }

      /**
       * @brief Returns the number of time intervals on which different sets of
       * basis vectors are defined.
       *
       * @return The number of time intervals on which there are basis vectors.
       */
      int
      getNumBasisTimeIntervals() const
      {
         return d_svd->getNumBasisTimeIntervals();
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
         CAROM_ASSERT(0 <= which_interval);
         CAROM_ASSERT(which_interval < getNumBasisTimeIntervals());
         return d_svd->getBasisIntervalStartTime(which_interval);
      }

      /**
       * @brief Returns true if the next state will result in a new time
       * interval.
       *
       * @return True if the next state results in the creation of a new time
       *         interval.
       */
      bool
      isNewTimeInterval() const
      {
         return d_svd->isNewTimeInterval();
      }

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      StaticSVDTimeStepper();

      /**
       * @brief Unimplemented copy constructor.
       */
      StaticSVDTimeStepper(
         const StaticSVDTimeStepper& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      StaticSVDTimeStepper&
      operator = (
         const StaticSVDTimeStepper& rhs);

      /**
       * @brief Pointer to the fundamental static SVD algorithm object.
       */
      boost::shared_ptr<StaticSVD> d_svd;
};

}

#endif
