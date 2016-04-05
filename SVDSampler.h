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

// Description: An abstract class defining the interface for determining the
//              next time at which a sample should be taken for basis
//              generation using an abstract SVD algorithm.

#ifndef included_SVDSampler_h
#define included_SVDSampler_h

#include "SVD.h"
#include <boost/shared_ptr.hpp>

namespace CAROM {

/**
 * Class SVDSampler defines the interface to determine the time at which the
 * next sample is needed when using an abstract SVD algorithm.
 */
class SVDSampler
{
   public:
      /**
       * @brief Constructor.
       */
      SVDSampler();

      /**
       * @brief Destructor.
       */
      ~SVDSampler();

      /**
       * @brief Returns true if it is time for the next sample.
       *
       * @param[in] time Time of interest.
       *
       * @return True if a sample should be taken at the supplied time.
       */
      virtual
      bool
      isNextSample(
         double time) = 0;

      /**
       * @brief Sample the new state, u_in.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] time The simulation time for the state.
       *
       * @return True if the sampling was successful.
       */
      bool
      takeSample(
         const double* u_in,
         double time)
      {
         return d_svd->takeSample(u_in, time);
      }

      /**
       * @brief Computes next time a state sample is needed.
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] rhs_in The right hand side at the specified time.
       * @param[in] time The simulation time for the state.
       *
       * @return The next time a statie sample is needed.
       */
      virtual
      double
      computeNextSampleTime(
         double* u_in,
         double* rhs_in,
         double time) = 0;

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

      /**
       * @brief Resets sample time step.
       *
       * @param[in] new_dt New value of sample time step.
       */
      virtual
      void
      resetDt(
         double new_dt) = 0;

      /**
       * @brief Returns the dimension of the system on this processor.
       *
       * @return The dimension of the system on this processor.
       */
      int
      getDim()
      {
         return d_svd->getDim();
      }

   protected:
      /**
       * @brief Pointer to the abstract SVD algorithm object.
       */
      boost::shared_ptr<SVD> d_svd;

   private:
      /**
       * @brief Unimplemented copy constructor.
       */
      SVDSampler(
         const SVDSampler& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      SVDSampler&
      operator = (
         const SVDSampler& rhs);
};

}

#endif
