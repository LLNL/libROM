/******************************************************************************
 *
 * Copyright (c) 2013-2016, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 * Written by William Arrighi wjarrighi@llnl.gov
 * CODE-686965
 * All rights reserved.
 *
 * This file is part of libROM.
 * For details, see https://computation.llnl.gov/librom
 * Please also read README_BSD_NOTICE.
 *
 * Redistribution and use in source and binary forms, with or without
 * modifications, are permitted provided that the following conditions are met:
 *
 *    o Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the disclaimer below.
 *    o Redistribution in binary form must reproduce the above copyright
 *      notice, this list of conditions and the disclaimer (as noted below) in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *    o Neither the name of the LLNS/LLNL nor the names of its contributors may
 *      be used to endorse or promote products derived from this software
 *      without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
 * LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OR SUCH DAMAGE.
 *
 *****************************************************************************/

// Description: The class that determines the next time at which a sample
//              should be taken for basis generation using the static SVD
//              approach.

#ifndef included_StaticSVDSampler_h
#define included_StaticSVDSampler_h

#include "SVDSampler.h"
#include "StaticSVD.h"
#include <boost/shared_ptr.hpp>

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
       * @param[in] debug_algorithm If true results of static svd algorithm
       *                            will be printed to facilitate debugging.
       */
      StaticSVDSampler(
         int dim,
         int samples_per_time_interval,
         bool debug_algorithm = false);

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
