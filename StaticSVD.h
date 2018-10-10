/******************************************************************************
 *
 * Copyright (c) 2013-2017, Lawrence Livermore National Security, LLC.
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

// Description: A class implementing interface of SVD for the static SVD
//              algorithm.

#ifndef included_StaticSVD_h
#define included_StaticSVD_h

#include "SVD.h"

namespace CAROM {

/**
 * StaticSVD implements the interface of class SVD for the static SVD
 * algorithm.  This algorithm is not scalable and is intended primarily as a
 * sanity check of the incremental svd algorithm.
 */
class StaticSVD : public SVD
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
      StaticSVD(
         int dim,
         int samples_per_time_interval,
         bool debug_algorithm = false);

      /**
       * Destructor.
       */
      ~StaticSVD();

      /**
       * @brief Collect the new sample, u_in at supplied time.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The new sample.
       * @param[in] time The simulation time of the new sample.
       * @param[in] add_without_increase If true, then addLinearlyDependent will be invoked
       *
       * @return True if the sampling was successful.
       */
      virtual
      bool
      takeSample(
         double* u_in,
         double time,
         bool add_without_increase = false);

      /**
       * @brief Returns the basis vectors for the current time interval.
       *
       * @post thisIntervalBasisCurrent()
       *
       * @return The basis vectors for the current time interval.
       */
      virtual
      const Matrix*
      getBasis();

      /**
       * @brief Returns the temporal basis vectors for the current time interval.
       *
       * @post thisIntervalBasisCurrent()
       *
       * @return The temporal basis vectors for the current time interval.
       */
      virtual
      const Matrix*
      getTBasis();

      /**
       * @brief Returns the singular values for the current time interval.
       *
       * @post thisIntervalBasisCurrent()
       *
       * @return The singular values for the current time interval.
       */
      virtual
      const Matrix*
      getSingularValues();

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      StaticSVD();

      /**
       * @brief Unimplemented copy constructor.
       */
      StaticSVD(
         const StaticSVD& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      StaticSVD&
      operator = (
         const StaticSVD& rhs);

      /**
       * @brief Gathers samples from all other processors to form complete
       * sample of system and computes the SVD.
       */
      void
      computeSVD();

      /**
       * @brief Preforms the actual SVD via lapack.
       *
       * Given a matrix, A, computes the 3 components of the singular value
       * decomposition.
       *
       * @pre A != 0
       * @pre total_dim > 0
       *
       * @param[in] A The globalized system whose SVD will be computed.
       * @param[in] total_dim The total dimension of the system that has been
       *                      distributed of multiple processors.
       */
      void
      svd(
         double* A,
         int total_dim);

      /**
       * @brief Tells if the basis vectors for this time interval are up to
       * date.
       *
       * @return True if the basis vectors for this time interval are up to
       * date.
       */
      bool
      thisIntervalBasisCurrent()
      {
         return d_this_interval_basis_current;
      }

      /**
       * @brief Current samples of the system.
       */
      std::vector<double*> d_samples;

      /**
       * @brief The globalized matrix L.
       *
       * L is small and each process owns all of L.
       */
      Matrix* d_V;

      /**
       * @brief Flag to indicate if the basis vectors for the current time
       * interval are up to date.
       */
      bool d_this_interval_basis_current;

      /**
       * @brief The rank of the process this object belongs to.
       */
      int d_rank;

      /**
       * @brief The number of processors being run on.
       */
      int d_num_procs;

      /**
       * @brief MPI message tag.
       */
      static const int COMMUNICATE_A;
};

}

#endif
