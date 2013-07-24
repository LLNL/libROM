#include "static_svd.h"
#include <boost/shared_ptr.hpp>

namespace CAROM {

// The class knows, given a static svd implementation, the time at which
// the next state collection is needed.  It also knows given a time whether
// it is time for the next state collection.
class static_svd_time_stepper
{
   public:
      // Constructor.
      static_svd_time_stepper(
         int* argc,
         char*** argv,
         int dim);

      // Destructor.
      ~static_svd_time_stepper();

      // Returns true if it is time for the next state collection.
      bool
      isNextStateCollection(
         double time)
      {
         return true;
      }

      // Collect the new state, u_in.
      void
      collectState(double* u_in)
      {
         d_svd->collectState(u_in);
      }

      // Computes next time a state collection is needed.
      double
      computeNextStateCollectionTime(
         double* u_in,
         double* rhs_in,
         double time)
      {
         return time;
      }

      // Returns the model parameters.
      Mat
      getModel()
      {
         return d_svd->getModel();
      }

   private:
      // Unimplemented default constructor.
      static_svd_time_stepper();

      // Unimplemented copy constructor.
      static_svd_time_stepper(const static_svd_time_stepper& other);

      // Unimplemented assignment operator.
      static_svd_time_stepper&
      operator = (const static_svd_time_stepper& rhs);

      // The fundamental static SVD algorithm.
      boost::shared_ptr<static_svd> d_svd;
};

}
