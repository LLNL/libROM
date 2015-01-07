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

#ifndef included_IncrementalSVDSampler_h
#define included_IncrementalSVDSampler_h

#include "SVDSampler.h"
#include "IncrementalSVD.h"
#include <boost/shared_ptr.hpp>

namespace CAROM {

/**
 * Class IncrementalSVDSampler knows, given an incremental svd implementation,
 * the time at which the next sample is needed.  It also knows given a time
 * whether it is time for the next sample.  There are two factors determining
 * if it is time for the next sample:
 * 1) the current time compared to the time the next sample must happen
 * 2) the number of time steps since the last sample
 */
class IncrementalSVDSampler : public SVDSampler
{
   public:
      /**
       * @brief Constructor.
       *
       * @pre dim > 0
       * @pre redundancy_tol > 0.0
       * @pre samples_per_time_interval > 0
       * @pre sampling_tol > 0.0
       * @pre max_time_between_samples > 0.0
       *
       * @param[in] dim The dimension of the system on this processor.
       * @param[in] redundancy_tol Tolerance to determine if a sample is
       *                           redundant or not.
       * @param[in] skip_redundant If true skip redundant samples.
       * @param[in] samples_per_time_interval The maximum number of samples in
       *                                      each time interval.
       * @param[in] sampling_tol Sampling control tolerance.  Limits error in
       *                         projection of solution into reduced order
       *                         space.
       * @param[in] max_time_between_samples Hard upper bound on time between
       *                                     samples.
       * @param[in] fast_update If true use the fast update incremental svd
       *                        algorithm.
       * @param[in] debug_algorithm If true results of incremental svd algorithm
       *                            will be printed to facilitate debugging.
       */
      IncrementalSVDSampler(
         int dim,
         double redundancy_tol,
         bool skip_redundant,
         int samples_per_time_interval,
         double sampling_tol,
         double max_time_between_samples,
         bool fast_update,
         bool debug_algorithm = false);

      /**
       * @brief Destructor.
       */
      ~IncrementalSVDSampler();

      /**
       * @brief Returns true if it is time for the next sample.
       *
       * @param[in] time Time of interest.
       *
       * @return True if it is time for the next sample to be taken.
       */
      virtual
      bool
      isNextSample(
         double time);

      /**
       * @brief Computes next time a sample is needed.
       *
       * @pre u_in != 0
       * @pre rhs_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] rhs_in The right hand side at the specified time.
       * @param[in] time The simulation time for the state.
       */
      virtual
      double
      computeNextSampleTime(
         double* u_in,
         double* rhs_in,
         double time);

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      IncrementalSVDSampler();

      /**
       * @brief Unimplemented copy constructor.
       */
      IncrementalSVDSampler(
         const IncrementalSVDSampler& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      IncrementalSVDSampler&
      operator = (
         const IncrementalSVDSampler& rhs);

      /**
       * @brief Sampling control tolerance.
       *
       * Limits error in projection of solution into the reduced order space.
       */
      double d_tol;

      /**
       * @brief Maximum time between samples.
       */
      double d_max_time_between_samples;

      /**
       * @brief Next time at which a sample should be taken.
       */
      double d_next_sample_time;
};

}

#endif
