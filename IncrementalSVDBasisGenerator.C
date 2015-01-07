/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
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
   double redundancy_tol,
   bool skip_redundant,
   int samples_per_time_interval,
   double sampling_tol,
   double max_time_between_samples,
   bool fast_update,
   const std::string& basis_file_name,
   bool debug_algorithm,
   Database::formats file_format) :
   SVDBasisGenerator(basis_file_name, file_format)
{
   CAROM_ASSERT(dim > 0);
   CAROM_ASSERT(redundancy_tol > 0.0);
   CAROM_ASSERT(samples_per_time_interval > 0);
   CAROM_ASSERT(sampling_tol > 0.0);
   CAROM_ASSERT(max_time_between_samples > 0.0);
   d_svdsampler.reset(new IncrementalSVDSampler(dim,
                                               redundancy_tol,
                                               skip_redundant,
                                               samples_per_time_interval,
                                               sampling_tol,
                                               max_time_between_samples,
                                               fast_update,
                                               debug_algorithm));
}

IncrementalSVDBasisGenerator::~IncrementalSVDBasisGenerator()
{
}

}
