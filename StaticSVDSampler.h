/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The class that determines the next time at which a sample
//              should be taken for basis generation using the static SVD
//              approach.

#ifndef included_StaticSVDSampler_h
#define included_StaticSVDSampler_h

#include "SVDSampler.h"
#include "StaticSVD.h"

namespace CAROM {

/**
 * Class StaticSVDSampler knows, given a static svd implementation, the
 * time at which the next sample is needed.  It also knows given a time whether
 * it is time for the next sample.  All state vectors are sampled in the
 * static SVD implementation so it is always time for a new sample.
 */
class StaticSVDSampler : public SVDSampler
{
   public:
      /**
       * @brief Constructor.
       *
       * @pre dim > 0
       * @pre samples_per_time_interval > 0
       *
       * @param[in] dim The dimension of the system on this processor.
       * @param[in] samples_per_time_interval The maximum number of samples
       *                                      in each time interval.
       * @param[in] max_time_intervals The maximium number of time intervals.
       * @param[in] max_basis_dimension (typemax(int)) The maximum number of
       *                                vectors returned in the basis.
       * @param[in] sigma_tolerance This tolerance is based on the ratio of
       *                            singular values to the largest singular
       *                            value. If sigma[i] / sigma[0] < sigma_tolerance,
       *                            the associated vector is dropped from the
       *                            basis.
       * @param[in] debug_algorithm If true results of static svd algorithm
       *                            will be printed to facilitate debugging.
       * @param[in] updateRightSV If true write right singular vectors
       */
      StaticSVDSampler(
         int dim,
         int samples_per_time_interval,
         int max_time_intervals = -1,
         int max_basis_dimension = std::numeric_limits<int>::max(),
         double sigma_tolerance = 0,
         bool debug_algorithm = false,
         bool updateRightSV = false);

      /**
       * @brief Destructor.
       */
      ~StaticSVDSampler();

      /**
       * @brief Returns true if it is time for the next sample.
       *
       * As the static algorithm samples everything this always returns true.
       *
       * @param[in] time Time of interest--unused.
       *
       * @return true
       */
      virtual
      bool
      isNextSample(
         double time);

      /**
       * @brief Computes next time a state sample is needed.
       *
       * @param[in] u_in The state at the specified time--unused.
       * @param[in] rhs_in The right hand side at the specified time--unused.
       * @param[in] time The simulation time for the state.
       *
       * @return The current simulation time as the static algorithm samples at
       * each time step.
       */
      virtual
      double
      computeNextSampleTime(
         double* u_in,
         double* rhs_in,
         double time);

      /**
       * @brief Resets sample time step.
       *
       * @param[in] new_dt New value of sample time step.
       */
      virtual
      void
      resetDt(
         double new_dt);

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      StaticSVDSampler();

      /**
       * @brief Unimplemented copy constructor.
       */
      StaticSVDSampler(
         const StaticSVDSampler& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      StaticSVDSampler&
      operator = (
         const StaticSVDSampler& rhs);
};

}

#endif
