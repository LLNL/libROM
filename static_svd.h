#include "petscmat.h"
#include "petscvec.h"
#include <vector>

namespace CAROM {

// A class which embodies the static SVD algorithm.  The API is intentionally
// small.  One may collect the states, compute the SVD, and get the dimension
// of the system.
class static_svd
{
   public:
      // Constructor.
      static_svd(
         int* argc,
         char*** argv,
         int dim);

      // Destructor.
      ~static_svd();

      // Collect the new state, u_in.
      void
      collectState(
         double* u_in)
      {
         d_state.push_back(u_in);
      }

      // Returns the dimension of the system.
      int
      getDim()
      {
         return d_dim;
      }

      // Returns the model parameters, a copy of d_U.
      Mat
      getModel()
      {
         computeSVD();
         Mat result;
         MatCreate(PETSC_COMM_WORLD, &result);
         MatSetSizes(result, d_dim, PETSC_DECIDE,
                     PETSC_DETERMINE, static_cast<int>(d_state.size()));
         MatSetType(result, MATDENSE);
         MatSetUp(result);
         MatCopy(d_U, result, SAME_NONZERO_PATTERN);
         return result;
      }

   private:
      // Unimplemented default constructor.
      static_svd();

      // Unimplemented copy constructor.
      static_svd(const static_svd& other);

      // Unimplemented assignment operator.
      static_svd&
      operator = (const static_svd& rhs);

      // Compute the SVD.
      void
      computeSVD();

      // Given a matrix, A, returns the 3 components of the singular value
      // decomposition.
      void
      svd(
         Mat A);

      // Dimension of the system.
      int d_dim;

      // Current state of the system.
      std::vector<double*> d_state;

      Mat d_U;
      Mat d_S;
      Mat d_V;

      // Rank of process this object lives on.
      int d_rank;
};

}
