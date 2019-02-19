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

// Description: The concrete implementation of the incremental SVD algorithm
//              that is equivalent to but computationally more expensive than
//              the "fast update" method.

#ifndef included_IncrementalSVDStandard_h
#define included_IncrementalSVDStandard_h

#include "IncrementalSVD.h"

namespace CAROM {

/**
 * A class which embodies the standard incremental SVD algorithm.
 */
class IncrementalSVDStandard : public IncrementalSVD
{
   public:
      /**
       * @brief Constructor.
       *
       * @pre dim > 0
       * @pre linearity_tol > 0.0
       * @pre samples_per_time_interval > 0
       *
       * @param[in] dim The dimension of the system on this processor.
       * @param[in] linearity_tol Tolerance to determine whether or not a
       *                          sample is linearly dependent.
       * @param[in] skip_linearly_dependent If true skip linearly dependent
       *                                    samples.
       * @param[in] samples_per_time_interval The number of samples to be
       *                                      collected for each time interval.
       * @param[in] basis_file_name The base part of the name of the file
       *                            containing the basis vectors.  Each process
       *                            will append its process ID to this base
       *                            name.
       * @param[in] save_state If true the state of the SVD will be written to
       *                       disk when the object is deleted.  If there are
       *                       multiple time intervals then the state will not
       *                       be saved as restoring such a state makes no
       *                       sense.
       * @param[in] restore_state If true the state of the SVD will be restored
       *                          when the object is created.
       * @param[in] debug_algorithm If true results of the algorithm will be
       *                            printed to facilitate debugging.
       */
      IncrementalSVDStandard(
         int dim,
         double linearity_tol,
         bool skip_linearly_dependent,
         int max_basis_dimension,
         int samples_per_time_interval,
         const std::string& basis_file_name,
         bool save_state = false,
         bool restore_state = false,
         bool updateRightSV = false,
         bool debug_algorithm = false);

      /**
       * @brief Destructor.
       */
      ~IncrementalSVDStandard();

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      IncrementalSVDStandard();

      /**
       * @brief Unimplemented copy constructor.
       */
      IncrementalSVDStandard(
         const IncrementalSVDStandard& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      IncrementalSVDStandard&
      operator = (
         const IncrementalSVDStandard& rhs);

      /**
       * @brief Constructs the first svd.
       *
       * @pre u != 0
       * @pre time >= 0.0
       *
       * @param[in] u The first state.
       * @param[in] time The simulation time for the first state.
       */
      virtual
      void
      buildInitialSVD(
         double* u,
         double time);

      /**
       * @brief Computes the current basis vectors.
       */
      virtual
      void
      computeBasis();

      /**
       * Add a linearly dependent sample to the svd.
       *
       * @pre A != 0
       * @pre sigma != 0
       *
       * @param[in] A The left singular vectors.
       * @param[in] W The right singular vectors.
       * @param[in] sigma The singular values.
       */
      void
      addLinearlyDependentSample(
         const Matrix* A,
         const Matrix* W,
         const Matrix* sigma);

      /**
       * @brief Add a new, unique sample to the svd.
       *
       * @pre j != 0
       * @pre A != 0
       * @pre W != 0
       * @pre sigma != 0
       *
       * @param[in] j The new column of d_U.
       * @param[in] A The left singular vectors.
       * @param[in] W The right singular vectors.
       * @param[in] sigma The singular values.
       */
      void
      addNewSample(
         const Vector* j,
         const Matrix* A,
         const Matrix* W,
         Matrix* sigma);
};

}

#endif
