#include "incremental_svd_time_stepper.h"

namespace CAROM {

// A class which implements a Reduced Order Model base on the incremental svd
// algorithm.  At the present time it pretty much just recasts the API of
// the incremental_svd_time_stepper in terms of snapshots.  In the future this
// class will likely need to do much more.
class incremental_svd_rom
{
   public:
      incremental_svd_rom(
         int* argc,
         char*** argv,
         int dim,
         double epsilon,
         bool skip_redundant,
         int max_time_steps_between_snapshots);

      ~incremental_svd_rom();

      bool
      isNextSnapshot(
         double time)
      {
         return d_isvdts->isNextIncrement(time);
      }

      void
      takeSnapshot(
         double* u_in)
      {
         d_isvdts->increment(u_in);
      }

      double
      computeNextSnapshotTime(
         double* u_in,
         double* rhs_in,
         double time)
      {
         return d_isvdts->computeNextIncrementTime(u_in, rhs_in, time);
      }

   private:
      boost::shared_ptr<incremental_svd_time_stepper> d_isvdts;
};

}
