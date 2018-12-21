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

// Description: The abstract wrapper class for an abstract SVD algorithm and
//              sampler.  This class provides interfaces to each so that an
//              application only needs to instantiate one concrete
//              implementation of this class to control all aspects of basis
//              vector generation.

#ifndef included_SVDBasisGenerator_h
#define included_SVDBasisGenerator_h

#include "BasisWriter.h"
#include "SVDSampler.h"

/* Use C++11 built-in shared pointers if available; else fallback to Boost. */
#if __cplusplus >= 201103L
#include <memory>
#else
#include <boost/shared_ptr.hpp>
#endif

#include <string.h>

namespace CAROM {

class BasisWriter;
class Matrix;

/**
 * Class SVDBasisGenerator is an abstract base class defining the interface for
 * the generation of basis vectors via the svd method.  This class wraps the
 * abstract SVD algorithm and sampler and provides interfaces to each so
 * that an application only needs to instantiate one concrete implementation of
 * this class to control all aspects of basis vector generation.
 */
class SVDBasisGenerator
{
   public:
      /**
       * @brief Destructor.
       */
      virtual
      ~SVDBasisGenerator();

      /**
       * @brief Returns true if it is time for the next svd sample.
       *
       * @pre time >= 0.0
       *
       * @param[in] time Time of interest.
       *
       * @return True if it is time for the next sample to be taken.
       */
      bool
      isNextSample(
         double time)
      {
         CAROM_ASSERT(time >= 0.0);
         return d_svdsampler->isNextSample(time);
      }
     
      /**
       * @brief Returns true if it needs to update right basis vectors.
       */
      bool
      updateRightSV() {return d_svdsampler->isUpdateRightSV(); };

      /**
       * @brief Sample the new state, u_in, at the given time.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] time The simulation time for the state.
       * @param[in] dt The current simulation dt.
       *
       * @return True if the sampling was successful.
       */
      bool
      takeSample(
         double* u_in,
         double time,
         double dt,
         bool add_without_increase = false)
      {
         CAROM_ASSERT(u_in != 0);
         CAROM_ASSERT(time >= 0);
 

         // Check that u_in is not non-zero.
         Vector u_vec(u_in, d_svdsampler->getDim(), true);
         if (u_vec.norm() == 0.0) {
            return false;
         }


         if (getNumBasisTimeIntervals() > 0 &&
             d_svdsampler->isNewTimeInterval()) {
            d_svdsampler->resetDt(dt);
            // YC: commenting this out unless someone think this is necessary
            //     we call writeBasis in "endSamples()" function below.
            //     I think that one call is enough. 
            //     I will remove completely if no one has other opinion.
            //if (d_basis_writer) {
            //   d_basis_writer->writeBasis();
            //}
         }

         return d_svdsampler->takeSample(u_in, time, add_without_increase);
      }

      /**
       * @brief Signal that the final sample has been taken.
       */
      void
      endSamples()
      {
         if (d_basis_writer) {
            d_basis_writer->writeBasis();
         }
      }

      /**
       * @brief Computes next time an svd sample is needed.
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
      computeNextSampleTime(
         double* u_in,
         double* rhs_in,
         double time)
      {
         CAROM_ASSERT(u_in != 0);
         CAROM_ASSERT(rhs_in != 0);
         CAROM_ASSERT(time >= 0);

         return d_svdsampler->computeNextSampleTime(u_in, rhs_in, time);
      }

      /**
       * @brief Returns the basis vectors for the current time interval as a
       * Matrix.
       *
       * @return The basis vectors for the current time interval.
       */
      const Matrix*
      getSpatialBasis()
      {
         return d_svdsampler->getSpatialBasis();
      }

      /**
       * @brief Returns the temporal basis vectors for the current time interval as a
       * Matrix.
       *
       * @return The temporal basis vectors for the current time interval.
       */
      const Matrix*
      getTemporalBasis()
      {
         return d_svdsampler->getTemporalBasis();
      }

      /**
       * @brief Returns the singular values for the current time interval as a
       * Matrix.
       *
       * @return The singular values for the current time interval.
       */
      const Matrix*
      getSingularValues()
      {
         return d_svdsampler->getSingularValues();
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
         return d_svdsampler->getNumBasisTimeIntervals();
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
         return d_svdsampler->getBasisIntervalStartTime(which_interval);
      }

   protected:
      /**
       * @brief Constructor.
       *
       * Although all member functions are implemented by delegation to either
       * d_basis_writer or d_svdsampler, this class is still abstract.  In this
       * context it is not yet known which SVDSampler to instantiate.  Hence an
       * instance of this class may not be constructed.
       *
       * @param[in] basis_file_name The base part of the name of the file
       *                            containing the basis vectors.  Each process
       *                            will append its process ID to this base
       *                            name.
       * @param[in] file_format The format of the file containing the basis
       *                        vectors.
       */
      SVDBasisGenerator(
         const std::string& basis_file_name,
         Database::formats file_format = Database::HDF5);

      /**
       * @brief Writer of basis vectors.
       */
      BasisWriter* d_basis_writer;

      /**
       * @brief Pointer to the underlying sampling control object.
       */
#if __cplusplus >= 201103L
      std::shared_ptr<SVDSampler> d_svdsampler;
#else
      boost::shared_ptr<SVDSampler> d_svdsampler;
#endif

   private:
      /**
       * @brief Unimplemented copy constructor.
       */
      SVDBasisGenerator(
         const SVDBasisGenerator& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      SVDBasisGenerator&
      operator = (
         const SVDBasisGenerator& rhs);
};

}

#endif
