/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The concrete wrapper class for static SVD algorithm and
//              sampler.  Implements interface of SVDBasisGenerator.

#include "StaticSVDBasisGenerator.h"
#include "StaticSVDSampler.h"

namespace CAROM {

StaticSVDBasisGenerator::StaticSVDBasisGenerator(
   StaticSVDBasisGeneratorOptions options,
   const std::string& basis_file_name,
   Database::formats file_format) :
   SVDBasisGenerator(basis_file_name, file_format)
{
   d_svdsampler.reset(new StaticSVDSampler(options.dim,
                                           options.samples_per_time_interval,
                                           options.max_basis_dimension,
                                           options.sigma_tolerance,
                                           options.debug_algorithm,
                                           options.output_rightSV));
}

StaticSVDBasisGenerator::~StaticSVDBasisGenerator()
{
}

}
