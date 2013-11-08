#ifndef included_incremental_svd_rom_h
#define included_incremental_svd_rom_h

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
         int dim,
         double epsilon,
         bool skip_redundant,
         int rom_size,
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

      // Add a snapshot to the incremental svd at the given time.
      void
      takeSnapshot(
         double* u_in,
         double time)
      {
         d_isvdts->increment(u_in, time);
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

      // Returns the model parameters at the given time as a Matrix.
      Matrix*
      getModel(
         double time)
      {
         return d_isvdts->getModel(time);
      }

      // Writes the model to a file with the given base name.
      void
      writeModel(
         const std::string& base_file_name)
      {
         d_isvdts->writeModel(base_file_name);
      }

      // Reads the model from a file with the given base name.
      void
      readModel(
         const std::string& base_file_name)
      {
         d_isvdts->readModel(base_file_name);
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

#endif
