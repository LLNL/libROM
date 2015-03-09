/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2015 Lawrence Livermore National Security, LLC
 * Description: The concrete wrapper class for a specific incremental SVD
 *              algorithm and sampler.  Implements interface of
 *              SVDBasisGenerator.
 *
 *****************************************************************************/

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
                                                min_sampling_time_step_scale,
                                                sampling_time_step_scale,
                                                max_sampling_time_step_scale,
                                                debug_algorithm));
}

IncrementalSVDBasisGenerator::~IncrementalSVDBasisGenerator()
{
}

}
