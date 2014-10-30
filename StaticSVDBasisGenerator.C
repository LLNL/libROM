/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The concrete wrapper class for static SVD algorithm and time
 *              stepper.  Implements interface of SVDBasisGenerator.
 *
 *****************************************************************************/

#include "StaticSVDBasisGenerator.h"

namespace CAROM {

StaticSVDBasisGenerator::StaticSVDBasisGenerator(
   int dim,
   int increments_per_time_interval,
   const std::string& basis_file_name,
   Database::formats file_format) :
   SVDBasisGenerator(basis_file_name, file_format),
   d_svdts(new StaticSVDTimeStepper(dim, increments_per_time_interval))
{
   CAROM_ASSERT(dim > 0);
   CAROM_ASSERT(increments_per_time_interval > 0);
}

StaticSVDBasisGenerator::~StaticSVDBasisGenerator()
{
}

bool
StaticSVDBasisGenerator::isNextSnapshot(
   double time)
{
   CAROM_ASSERT(time >= 0.0);
   return d_svdts->isNextStateCollection(time);
}

void
StaticSVDBasisGenerator::takeSnapshot(
   double* u_in,
   double time)
{
   CAROM_ASSERT(u_in != 0);
   CAROM_ASSERT(time >= 0);

   if (d_basis_writer &&
       d_svdts->isNewTimeInterval() && getNumBasisTimeIntervals() > 0) {
      d_basis_writer->writeBasis();
   }
   d_svdts->collectState(u_in, time);
}

double
StaticSVDBasisGenerator::computeNextSnapshotTime(
   double* u_in,
   double* rhs_in,
   double time)
{
   CAROM_ASSERT(u_in != 0);
   CAROM_ASSERT(rhs_in != 0);
   CAROM_ASSERT(time >= 0);

   return d_svdts->computeNextStateCollectionTime(u_in, rhs_in, time);
}

const Matrix*
StaticSVDBasisGenerator::getBasis()
{
   return d_svdts->getBasis();
}

int
StaticSVDBasisGenerator::getNumBasisTimeIntervals() const
{
   return d_svdts->getNumBasisTimeIntervals();
}

double
StaticSVDBasisGenerator::getBasisIntervalStartTime(
   int which_interval) const
{
   CAROM_ASSERT(0 <= which_interval);
   CAROM_ASSERT(which_interval < getNumBasisTimeIntervals());
   return d_svdts->getBasisIntervalStartTime(which_interval);
}

}
