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
   int basis_size,
   double sampling_tol,
   double max_time_between_snapshots,
   bool fast_update,
   const std::string& basis_file_name,
   Database::formats file_format) :
   SVDBasisGenerator(basis_file_name, file_format),
   d_isvdts(new IncrementalSVDTimeStepper(dim,
                                          redundancy_tol,
                                          skip_redundant,
                                          basis_size,
                                          sampling_tol,
                                          max_time_between_snapshots,
                                          fast_update))
{
}

IncrementalSVDBasisGenerator::~IncrementalSVDBasisGenerator()
{
}

bool
IncrementalSVDBasisGenerator::isNextSnapshot(
   double time)
{
   return d_isvdts->isNextIncrement(time);
}

void
IncrementalSVDBasisGenerator::takeSnapshot(
   double* u_in,
   double time)
{
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
   return d_isvdts->getBasisIntervalStartTime(which_interval);
}

}
