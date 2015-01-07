/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The concrete wrapper class for static SVD algorithm and
 *              sampler.  Implements interface of SVDBasisGenerator.
 *
 *****************************************************************************/

#include "StaticSVDBasisGenerator.h"
#include "StaticSVDSampler.h"

namespace CAROM {

StaticSVDBasisGenerator::StaticSVDBasisGenerator(
   int dim,
   int samples_per_time_interval,
   const std::string& basis_file_name,
   bool debug_algorithm,
   Database::formats file_format) :
   SVDBasisGenerator(basis_file_name, file_format)
{
   CAROM_ASSERT(dim > 0);
   CAROM_ASSERT(samples_per_time_interval > 0);
   d_svdsampler.reset(new StaticSVDSampler(dim,
                                           samples_per_time_interval,
                                           debug_algorithm));
}

StaticSVDBasisGenerator::~StaticSVDBasisGenerator()
{
}

}
