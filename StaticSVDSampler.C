/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The class that determines the next time at which a sample
 *              should be taken for basis generation using the static SVD
 *              approach.
 *
 *****************************************************************************/

#include "StaticSVDSampler.h"

namespace CAROM {

StaticSVDSampler::StaticSVDSampler(
   int dim,
   int samples_per_time_interval) :
   d_svd(new StaticSVD(dim, samples_per_time_interval))
{
   CAROM_ASSERT(dim > 0);
   CAROM_ASSERT(samples_per_time_interval > 0);
}

StaticSVDSampler::~StaticSVDSampler()
{
}

}
