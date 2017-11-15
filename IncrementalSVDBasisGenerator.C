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

#include "IncrementalSVDBasisGenerator.h"
#include "IncrementalSVDSampler.h"

namespace CAROM {

IncrementalSVDBasisGenerator::IncrementalSVDBasisGenerator(
   int dim,
   double linearity_tol,
   bool skip_linearly_dependent,
   bool fast_update,
   double initial_dt,
   int samples_per_time_interval,
   double sampling_tol,
   double max_time_between_samples,
   const std::string& basis_file_name,
   bool save_state,
   bool restore_state,
   Database::formats file_format,
   double min_sampling_time_step_scale,
   double sampling_time_step_scale,
   double max_sampling_time_step_scale,
   bool debug_algorithm) :
   SVDBasisGenerator(basis_file_name, file_format)
{
   CAROM_ASSERT(dim > 0);
   CAROM_ASSERT(linearity_tol > 0.0);
   CAROM_ASSERT(initial_dt > 0.0);
   CAROM_ASSERT(samples_per_time_interval > 0);
   CAROM_ASSERT(sampling_tol > 0.0);
   CAROM_ASSERT(max_time_between_samples > 0.0);
   CAROM_ASSERT(min_sampling_time_step_scale < max_sampling_time_step_scale);

   d_svdsampler.reset(new IncrementalSVDSampler(dim,
                                                linearity_tol,
                                                skip_linearly_dependent,
                                                fast_update,
                                                initial_dt,
                                                samples_per_time_interval,
                                                sampling_tol,
                                                max_time_between_samples,
                                                save_state,
                                                restore_state,
                                                min_sampling_time_step_scale,
                                                sampling_time_step_scale,
                                                max_sampling_time_step_scale,
                                                debug_algorithm));
}

IncrementalSVDBasisGenerator::~IncrementalSVDBasisGenerator()
{
}

}
