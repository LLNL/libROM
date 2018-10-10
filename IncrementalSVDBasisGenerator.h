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

// Description: The concrete wrapper class for a specific incremental SVD
//              algorithm and sampler.  Implements interface of
//              SVDBasisGenerator.

#ifndef included_IncrementalSVDBasisGenerator_h
#define included_IncrementalSVDBasisGenerator_h

#include "SVDBasisGenerator.h"

namespace CAROM {

/**
 * Class IncrementalSVDBasisGenerator implements the interface of base class
 * SVDBasisGenerator for the incremental svd algorithm.  Either the fast update
 * or the standard incremental algorithm may be specified through the
 * constructor.
 */
class IncrementalSVDBasisGenerator : public SVDBasisGenerator
{
   public:
      /**
       * @brief Constructor.
       *
       * @pre dim > 0
       * @pre linearity_tol > 0.0
       * @pre initial_dt > 0.0
       * @pre samples_per_time_interval > 0
       * @pre sampling_tol > 0.0
       * @pre max_time_between_samples > 0.0
       * @pre min_sampling_time_step_scale >= 0.0
       * @pre sampling_time_step_scale >= 0.0
       * @pre max_sampling_time_step_scale >= 0.0
       * @pre min_sampling_time_step_scale <= max_sampling_time_step_scale
       *
       * @param[in] dim The dimension of the system on this processor.
       * @param[in] linearity_tol Tolerance to determine whether or not a
       *                          sample is linearly dependent.
       * @param[in] skip_linearly_dependent If true skip linearly dependent
       *                                    samples.
       * @param[in] fast_update If true use the fast update algorithm.
       * @param[in] initial_dt Initial simulation time step.
       * @param[in] samples_per_time_interval The maximum number of samples in
       *                                      each time interval.
       * @param[in] sampling_tol Sampling control tolerance.  Limits error in
       *                         projection of solution into reduced order
       *                         space.
       * @param[in] max_time_between_samples Upper bound on time between
       *                                     samples.
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
       * @param[in] file_format The format of the file containing the basis
       *                        vectors.
       * @param[in] min_sampling_time_step_scale Minimum overall scale factor
       *                                         to apply to time step.
       * @param[in] sampling_time_step_scale Scale factor to apply to sampling
       *                                     algorithm.
       * @param[in] max_sampling_time_step_scale Maximum overall scale factor
       *                                         to apply to time step.
       * @param[in] debug_algorithm If true results of incremental svd
       *                            algorithm will be printed to facilitate
       *                            debugging.
       */
      IncrementalSVDBasisGenerator(
         int dim,
         double linearity_tol,
         bool skip_linearly_dependent,
         bool fast_update,
         int max_basis_dimension,
         double initial_dt,
         int samples_per_time_interval,
         double sampling_tol,
         double max_time_between_samples,
         const std::string& basis_file_name = "",
         bool save_state = false,
         bool restore_state = false,
         bool updateRightSV = false, 
         Database::formats file_format = Database::HDF5,
         double min_sampling_time_step_scale = 0.1,
         double sampling_time_step_scale = 0.8,
         double max_sampling_time_step_scale = 5.0,
         bool debug_algorithm = false);

      /**
       * @brief Destructor.
       */
      virtual
      ~IncrementalSVDBasisGenerator();

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      IncrementalSVDBasisGenerator();

      /**
       * @brief Unimplemented copy constructor.
       */
      IncrementalSVDBasisGenerator(
         const IncrementalSVDBasisGenerator& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      IncrementalSVDBasisGenerator&
      operator = (
         const IncrementalSVDBasisGenerator& rhs);
};

}

#endif
