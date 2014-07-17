#include "incremental_svd_basis_generator.h"

namespace CAROM {

incremental_svd_basis_generator::incremental_svd_basis_generator(
   int dim,
   double redundancy_tol,
   bool skip_redundant,
   int basis_size,
   double sampling_tol,
   double max_time_between_snapshots,
   bool fast_update,
   const std::string& basis_file_name,
   Database::formats file_format) :
   svd_basis_generator(basis_file_name, file_format),
   d_isvdts(new incremental_svd_time_stepper(dim,
                                             redundancy_tol,
                                             skip_redundant,
                                             basis_size,
                                             sampling_tol,
                                             max_time_between_snapshots,
                                             fast_update))
{
}

incremental_svd_basis_generator::~incremental_svd_basis_generator()
{
}

bool
incremental_svd_basis_generator::isNextSnapshot(
   double time)
{
   return d_isvdts->isNextIncrement(time);
}

void
incremental_svd_basis_generator::takeSnapshot(
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
incremental_svd_basis_generator::computeNextSnapshotTime(
   double* u_in,
   double* rhs_in,
   double time)
{
   return d_isvdts->computeNextIncrementTime(u_in, rhs_in, time);
}

const Matrix*
incremental_svd_basis_generator::getBasis()
{
   return d_isvdts->getBasis();
}

int
incremental_svd_basis_generator::getNumBasisTimeIntervals() const
{
   return d_isvdts->getNumBasisTimeIntervals();
}

double
incremental_svd_basis_generator::getBasisIntervalStartTime(
   int which_interval) const
{
   return d_isvdts->getBasisIntervalStartTime(which_interval);
}

}
