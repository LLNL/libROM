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

namespace CAROM {

StaticSVDBasisGenerator::StaticSVDBasisGenerator(
   int dim,
   int samples_per_time_interval,
   const std::string& basis_file_name,
   bool debug_rom,
   Database::formats file_format) :
   SVDBasisGenerator(basis_file_name, file_format),
   d_svdsampler(new StaticSVDSampler(dim,
                                     samples_per_time_interval,
                                     debug_rom))
{
   CAROM_ASSERT(dim > 0);
   CAROM_ASSERT(samples_per_time_interval > 0);
}

StaticSVDBasisGenerator::~StaticSVDBasisGenerator()
{
}

bool
StaticSVDBasisGenerator::isNextSample(
   double time)
{
   CAROM_ASSERT(time >= 0.0);
   return d_svdsampler->isNextSample(time);
}

void
StaticSVDBasisGenerator::takeSample(
   double* u_in,
   double time)
{
   CAROM_ASSERT(u_in != 0);
   CAROM_ASSERT(time >= 0);

   if (d_basis_writer &&
       d_svdsampler->isNewTimeInterval() && getNumBasisTimeIntervals() > 0) {
      d_basis_writer->writeBasis();
   }
   d_svdsampler->takeSample(u_in, time);
}

double
StaticSVDBasisGenerator::computeNextSampleTime(
   double* u_in,
   double* rhs_in,
   double time)
{
   CAROM_ASSERT(u_in != 0);
   CAROM_ASSERT(rhs_in != 0);
   CAROM_ASSERT(time >= 0);

   return d_svdsampler->computeNextSampleTime(u_in, rhs_in, time);
}

const Matrix*
StaticSVDBasisGenerator::getBasis()
{
   return d_svdsampler->getBasis();
}

int
StaticSVDBasisGenerator::getNumBasisTimeIntervals() const
{
   return d_svdsampler->getNumBasisTimeIntervals();
}

double
StaticSVDBasisGenerator::getBasisIntervalStartTime(
   int which_interval) const
{
   CAROM_ASSERT(0 <= which_interval);
   CAROM_ASSERT(which_interval < getNumBasisTimeIntervals());
   return d_svdsampler->getBasisIntervalStartTime(which_interval);
}

}
