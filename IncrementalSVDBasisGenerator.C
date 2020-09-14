/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The concrete wrapper class for a specific incremental SVD
//              algorithm and sampler.  Implements interface of
//              SVDBasisGenerator.

#include "IncrementalSVDBasisGenerator.h"
#include "IncrementalSVDSampler.h"

namespace CAROM {

IncrementalSVDBasisGenerator::IncrementalSVDBasisGenerator(
   IncrementalSVDBasisGeneratorOptions options,
   const std::string& basis_file_name,
   Database::formats file_format) :
   SVDBasisGenerator(basis_file_name, file_format)
{
   d_svdsampler.reset(new IncrementalSVDSampler(options.dim,
                                                options.linearity_tol,
                                                options.skip_linearly_dependent,
                                                options.fast_update,
                                                options.max_basis_dimension,
                                                options.initial_dt,
                                                options.samples_per_time_interval,
                                                options.sampling_tol,
                                                options.max_time_between_samples,
                                                basis_file_name,
                                                options.save_state,
                                                options.restore_state,
                                                options.updateRightSV,
                                                options.min_sampling_time_step_scale,
                                                options.sampling_time_step_scale,
                                                options.max_sampling_time_step_scale,
                                                options.debug_algorithm));
}

IncrementalSVDBasisGenerator::~IncrementalSVDBasisGenerator()
{
}

}
