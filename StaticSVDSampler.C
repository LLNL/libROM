/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The class that determines the next time at which a sample
//              should be taken for basis generation using the static SVD
//              approach.

#include "StaticSVDSampler.h"

namespace CAROM {

StaticSVDSampler::StaticSVDSampler(
   StaticSVDOptions options)
{
   d_svd.reset(new StaticSVD(options));
   d_updateRightSV = options.output_rightSV;
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

void
StaticSVDSampler::resetDt(
   double new_dt)
{
   CAROM_NULL_USE(new_dt);
}

}
