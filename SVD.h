/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC.
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

// Description: An abstract class defining the interface to the generic SVD
//              algorithm.

#ifndef included_SVD_h
#define included_SVD_h

#include "Matrix.h"
#include <vector>

namespace CAROM {

/**
 * Class SVD defines the interface to the generic SVD algorithm.  The API is
 * intentionally small.  One may collect the samples, compute the SVD, and get
 * the dimension of the system.
 */
class SVD
{
   public:
      /**
       * @brief Constructor.
       *
       * @pre dim > 0
       * @pre sample_per_time_interval > 0
       *
       * @param[in] dim The dimension of the system distributed to this
       *                processor.
       * @param[in] samples_per_time_interval The maximium number of samples
       *                                      collected in a time interval.
       * @param[in] debug_algorithm If true results of the algorithm will be
       *                            printed to facilitate debugging.
       */
      SVD(
         int dim,
         int samples_per_time_interval,
         bool debug_algorithm = false);

      /**
       * Destructor.
       */
      ~SVD();

      /**
       * @brief Collect the new sample, u_in at supplied time.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The new sample.
       * @param[in] time The simulation time of the new sample.
       * @param[in] add_without_increase If true, the addLinearlyDependent is invoked. 
       *
       * @return True if the sampling was successful.
       */
      virtual
      bool
      takeSample(
         double* u_in,
         double time,
         bool add_without_increase) = 0;

      /**
       * @brief Returns the dimension of the system on this processor.
       *
       * @return The dimension of the system on this processor.
       */
      int
      getDim() const
      {
         return d_dim;
      }

      /**
       * @brief Returns the basis vectors for the current time interval.
       *
       * @return The basis vectors for the current time interval.
       */
      virtual
      const Matrix*
      getSpatialBasis() = 0;

      /**
       * @brief Returns the temporal basis vectors for the current time interval.
       *
       * @return The temporal basis vectors for the current time interval.
       */
      virtual
      const Matrix*
      getTemporalBasis() = 0;

      /**
       * @brief Returns the singular values for the current time interval.
       *
       * @return The singular values for the current time interval.
       */
      virtual
      const Matrix*
      getSingularValues() = 0;

      /**
       * @brief Returns the number of time intervals on which different sets
       * of basis vectors are defined.
       *
       * @return The number of time intervals on which there are basis vectors.
       */
      int
      getNumBasisTimeIntervals() const
      {
         return static_cast<int>(d_time_interval_start_times.size());
      }

      /**
       * @brief Returns the start time for the requested time interval.
       *
       * @pre 0 <= which_interval
       * @pre which_interval < getNumBasisTimeIntervals()
       *
       * @param[in] which_interval The time interval of interest.
       *
       * @return The start time for the requested time interval.
       */
      double
      getBasisIntervalStartTime(
         int which_interval) const
      {
         CAROM_ASSERT(0 <= which_interval);
         CAROM_ASSERT(which_interval < getNumBasisTimeIntervals());

	 std::size_t i = static_cast<std::size_t>(which_interval);
         return d_time_interval_start_times[i];
      }

      /**
       * @brief Returns true if the next sample will result in a new time
       * interval.
       *
       * @return True if the next sample results in the creation of a new time
       *         interval.
       */
      bool
      isNewTimeInterval() const
      {
         return (d_num_samples == 0) ||
                (d_num_samples >= d_samples_per_time_interval);
      }

   protected:
      /**
       * @brief Dimension of the system.
       */
      int d_dim;

      /**
       * @brief Number of samples stored for the current time interval.
       */
      int d_num_samples;

      /**
       * @brief Number of rows in right singular matrix.
       */
      int d_num_rows_of_W;

      /**
       * @brief The maximum number of samples to be collected for a time
       * interval.
       */
      int d_samples_per_time_interval;

      /**
       * @brief The globalized basis vectors for the current time interval.
       *
       * The basis vectors are large and each process owns all of the basis
       * vectors.
       */
      Matrix* d_basis;

      /**
       * @brief The globalized right basis vectors for the current time interval.
       *
       * Depending on the SVD algorithm, it may be  distributed across all
       * processors or each processor may own all of U.
       */
      Matrix* d_basis_right;

      /**
       * @brief The matrix U which is large.
       *
       * Depending on the SVD algorithm, d_U may be  distributed across all
       * processors or each processor may own all of U.
       */
      Matrix* d_U;

      /**
       * @brief The matrix U which is large.
       *
       * Depending on the SVD algorithm, d_W may be  distributed across all
       * processors or each processor may own all of U.
       */
      Matrix* d_W;

      /**
       * @brief The matrix S which is small.
       *
       * For all SVD algorithms, S is not distributed and the entire matrix
       * exists on each processor.
       */
      Matrix* d_S;

      /**
       * @brief The simulation time at which each time interval starts.
       */
      std::vector<double> d_time_interval_start_times;

      /**
       * @brief Flag to indicate if results of algorithm should be printed for
       * debugging purposes.
       */
      bool d_debug_algorithm;

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      SVD();

      /**
       * @brief Unimplemented copy constructor.
       */
      SVD(
         const SVD& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      SVD&
      operator = (
         const SVD& rhs);
};

}

#endif
