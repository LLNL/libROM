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
   int samples_per_time_interval,
   bool debug_rom)
{
   CAROM_ASSERT(dim > 0);
   CAROM_ASSERT(samples_per_time_interval > 0);
   d_svd.reset(new StaticSVD(dim, samples_per_time_interval, debug_rom));
}

StaticSVDSampler::~StaticSVDSampler()
{
}

bool
StaticSVDSampler::isNextSample(
   double time)
{
   CAROM_NULL_USE(time);
   return true;
}

double
StaticSVDSampler::computeNextSampleTime(
   double* u_in,
   double* rhs_in,
   double time)
{
   CAROM_NULL_USE(u_in);
   CAROM_NULL_USE(rhs_in);
   return time;
}

}
