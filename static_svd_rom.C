#include "static_svd_rom.h"

namespace CAROM {

static_svd_rom::static_svd_rom(
   int dim,
   int increments_per_time_interval,
   const std::string& basis_file_name,
   Database::formats file_format) :
   svd_rom(basis_file_name, file_format),
   d_svdts(new static_svd_time_stepper(dim, increments_per_time_interval))
{
}

static_svd_rom::~static_svd_rom()
{
}

bool
static_svd_rom::isNextSnapshot(
   double time)
{
   return d_svdts->isNextStateCollection(time);
}

void
static_svd_rom::takeSnapshot(
   double* u_in,
   double time)
{
   if (d_basis_writer &&
       d_svdts->isNewTimeInterval() && getNumBasisTimeIntervals() > 0) {
      d_basis_writer->writeBasis();
   }
   d_svdts->collectState(u_in, time);
}

double
static_svd_rom::computeNextSnapshotTime(
   double* u_in,
   double* rhs_in,
   double time)
{
   return d_svdts->computeNextStateCollectionTime(u_in, rhs_in, time);
}

const Matrix*
static_svd_rom::getBasis()
{
   return d_svdts->getBasis();
}

int
static_svd_rom::getNumBasisTimeIntervals() const
{
   return d_svdts->getNumBasisTimeIntervals();
}

double
static_svd_rom::getBasisIntervalStartTime(
   int which_interval) const
{
   return d_svdts->getBasisIntervalStartTime(which_interval);
}

}
