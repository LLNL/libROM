#include "incremental_svd_rom.h"

namespace CAROM {

incremental_svd_rom::incremental_svd_rom(
   int dim,
   double epsilon,
   bool skip_redundant,
   int rom_size,
   double tolerance,
   double max_time_between_snapshots,
   bool fast_update,
   const std::string& basis_file_name,
   Database::formats file_format) :
   svd_rom(basis_file_name, file_format),
   d_isvdts(new incremental_svd_time_stepper(dim,
                                             epsilon,
                                             skip_redundant,
                                             rom_size,
                                             tolerance,
                                             max_time_between_snapshots,
                                             fast_update))
{
}

incremental_svd_rom::~incremental_svd_rom()
{
}

bool
incremental_svd_rom::isNextSnapshot(
   double time)
{
   return d_isvdts->isNextIncrement(time);
}

void
incremental_svd_rom::takeSnapshot(
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
incremental_svd_rom::computeNextSnapshotTime(
   double* u_in,
   double* rhs_in,
   double time)
{
   return d_isvdts->computeNextIncrementTime(u_in, rhs_in, time);
}

const Matrix*
incremental_svd_rom::getBasis()
{
   return d_isvdts->getBasis();
}

int
incremental_svd_rom::getNumBasisTimeIntervals() const
{
   return d_isvdts->getNumBasisTimeIntervals();
}

double
incremental_svd_rom::getBasisIntervalStartTime(
   int which_interval) const
{
   return d_isvdts->getBasisIntervalStartTime(which_interval);
}

}
