#include "incremental_svd_time_stepper.h"

namespace CAROM {

// A class which implements a Reduced Order Model base on the incremental svd
// algorithm.  At the present time it pretty much just recasts the API of
// the incremental_svd_time_stepper in terms of snapshots.  In the future this
// class will likely need to do much more.
class incremental_svd_rom
{
   public:
      // Constructor.
      incremental_svd_rom(
         int* argc,
         char*** argv,
         int dim,
         double epsilon,
         bool skip_redundant,
         int max_time_steps_between_snapshots);

      // Destructor.
      ~incremental_svd_rom();

      // Returns true if it is time for the next svd snapshot.
      bool
      isNextSnapshot(
         double time)
      {
         return d_isvdts->isNextIncrement(time);
      }

      // Add a snapshot to the incremental svd.
      void
      takeSnapshot(
         double* u_in)
      {
         d_isvdts->increment(u_in);
      }

      // Computes next time an svd snapshot is needed.
      double
      computeNextSnapshotTime(
         double* u_in,
         double* rhs_in,
         double time)
      {
         return d_isvdts->computeNextIncrementTime(u_in, rhs_in, time);
      }

      // Returns the model parameters as a PETSc Mat.
      Mat
      getModel()
      {
         return d_isvdts->getModel();
      }

   private:
      // Unimplemented default constructor.
      incremental_svd_rom();

      // Unimplemented copy constructor.
      incremental_svd_rom(const incremental_svd_rom& other);

      // Unimplemented assignment operator.
      incremental_svd_rom&
      operator = (const incremental_svd_rom& rhs);

      // The underlying time step control object.
      boost::shared_ptr<incremental_svd_time_stepper> d_isvdts;
};

}
