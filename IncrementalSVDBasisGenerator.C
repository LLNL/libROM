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
   int dim,
   double linearity_tol,
   bool skip_linearly_dependent,
   bool fast_update,
   int max_basis_dimension,
   double initial_dt,
   int samples_per_time_interval,
   double sampling_tol,
   double max_time_between_samples,
   int max_time_intervals,
   const std::string& basis_file_name,
   bool save_state,
   bool restore_state,
   bool updateRightSV,
   bool write_snapshots,
   Database::formats file_format,
   double min_sampling_time_step_scale,
   double sampling_time_step_scale,
   double max_sampling_time_step_scale,
   bool debug_algorithm) :
   SVDBasisGenerator(basis_file_name, file_format, write_snapshots)
{
   d_svdsampler.reset(new IncrementalSVDSampler(dim,
                                                linearity_tol,
                                                skip_linearly_dependent,
                                                fast_update,
                                                max_basis_dimension,
                                                initial_dt,
                                                samples_per_time_interval,
                                                sampling_tol,
                                                max_time_between_samples,
                                                max_time_intervals,
                                                basis_file_name,
                                                save_state,
                                                restore_state,
                                                updateRightSV,
                                                min_sampling_time_step_scale,
                                                sampling_time_step_scale,
                                                max_sampling_time_step_scale,
                                                debug_algorithm));
}

IncrementalSVDBasisGenerator::~IncrementalSVDBasisGenerator()
{
}

}
