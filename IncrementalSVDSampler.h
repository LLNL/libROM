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
//              should be taken for basis generation using an incremental SVD
//              approach.

#ifndef included_IncrementalSVDSampler_h
#define included_IncrementalSVDSampler_h

#include "SVDSampler.h"
#include "IncrementalSVD.h"

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
       * @param[in] options The struct containing the options for this basis
       *                    generator.
       * @param[in] basis_file_name The base part of the name of the file
       *                            containing the basis vectors.  Each process
       *                            will append its process ID to this base
       *                            name.
       */
      IncrementalSVDSampler(
         IncrementalSVDOptions options,
         const std::string& basis_file_name = "");

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
       * @brief Minimum sampling time step scale factor.
       */
      double d_min_sampling_time_step_scale;

      /**
       * @brief Sampling time step scale factor to apply to algorithm.
       */
      double d_sampling_time_step_scale;

      /**
       * @brief Maximum sampling time step scale factor.
       */
      double d_max_sampling_time_step_scale;

      /**
       * @brief Current time step.
       */
      double d_dt;

      /**
       * @brief Next time at which a sample should be taken.
       */
      double d_next_sample_time;

      /**
       * @brief The number of processors being run on.
       */
      int d_num_procs;
};

}

#endif
