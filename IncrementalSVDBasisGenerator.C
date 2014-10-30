/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The concrete wrapper class for a specific incremental SVD
 *              algorithm and time stepper.  Implements interface of
 *              SVDBasisGenerator.
 *
 *****************************************************************************/

#include "IncrementalSVDBasisGenerator.h"

namespace CAROM {

IncrementalSVDBasisGenerator::IncrementalSVDBasisGenerator(
   int dim,
   double redundancy_tol,
   bool skip_redundant,
   int increments_per_time_interval,
   double sampling_tol,
   double max_time_between_snapshots,
   bool fast_update,
   const std::string& basis_file_name,
   Database::formats file_format) :
   SVDBasisGenerator(basis_file_name, file_format),
   d_isvdts(new IncrementalSVDTimeStepper(dim,
                                          redundancy_tol,
                                          skip_redundant,
                                          increments_per_time_interval,
                                          sampling_tol,
                                          max_time_between_snapshots,
                                          fast_update))
{
   CAROM_ASSERT(dim > 0);
   CAROM_ASSERT(redundancy_tol > 0.0);
   CAROM_ASSERT(increments_per_time_interval > 0);
   CAROM_ASSERT(sampling_tol > 0.0);
   CAROM_ASSERT(max_time_between_snapshots > 0.0);
}

IncrementalSVDBasisGenerator::~IncrementalSVDBasisGenerator()
{
}

bool
IncrementalSVDBasisGenerator::isNextSnapshot(
   double time)
{
   CAROM_ASSERT(time >= 0);
   return d_isvdts->isNextIncrement(time);
}

void
IncrementalSVDBasisGenerator::takeSnapshot(
   double* u_in,
   double time)
{
   CAROM_ASSERT(u_in != 0);
   CAROM_ASSERT(time >= 0);

   if (d_basis_writer &&
       d_isvdts->isNewTimeInterval() && getNumBasisTimeIntervals() > 0) {
      d_basis_writer->writeBasis();
   }
   d_isvdts->increment(u_in, time);
}

double
IncrementalSVDBasisGenerator::computeNextSnapshotTime(
   double* u_in,
   double* rhs_in,
   double time)
{
   CAROM_ASSERT(u_in != 0);
   CAROM_ASSERT(rhs_in != 0);
   CAROM_ASSERT(time >= 0);

   return d_isvdts->computeNextIncrementTime(u_in, rhs_in, time);
}

const Matrix*
IncrementalSVDBasisGenerator::getBasis()
{
   return d_isvdts->getBasis();
}

int
IncrementalSVDBasisGenerator::getNumBasisTimeIntervals() const
{
   return d_isvdts->getNumBasisTimeIntervals();
}

double
IncrementalSVDBasisGenerator::getBasisIntervalStartTime(
   int which_interval) const
{
   CAROM_ASSERT(0 <= which_interval);
   CAROM_ASSERT(which_interval < getNumBasisTimeIntervals());
   return d_isvdts->getBasisIntervalStartTime(which_interval);
}

}
