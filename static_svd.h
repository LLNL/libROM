#ifndef included_static_svd_h
#define included_static_svd_h

#include "matrix.h"
#include <vector>

namespace CAROM {

// A class which embodies the static SVD algorithm.  The API is intentionally
// small.  One may collect the states, compute the SVD, and get the dimension
// of the system.
// This implementation is not scalable and is intended primarily as a sanity
// check of the incremental svd algorithm.
class static_svd
{
   public:
      // Constructor.
      static_svd(
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
      getDim() const
      {
         return d_dim;
      }

      // Returns the model parameters, d_U.
      const Matrix*
      getModel()
      {
         computeSVD();
         return d_U;
      }

   private:
      // Unimplemented default constructor.
      static_svd();

      // Unimplemented copy constructor.
      static_svd(
         const static_svd& other);

      // Unimplemented assignment operator.
      static_svd&
      operator = (
         const static_svd& rhs);

      // Compute the SVD.
      void
      computeSVD();

      // Given a matrix, A, returns the 3 components of the singular value
      // decomposition.
      void
      svd(
         double* A);

      // Dimension of the system.
      int d_dim;

      // Current state of the system.
      std::vector<double*> d_state;

      // The globalized matrix U.  U is large and each process owns all of U.
      Matrix* d_U;

      // The globalized matrix S.  S is small and each process owns all of S.
      Matrix* d_S;

      // The globalized matrix L.  L is small and each process owns all of L.
      Matrix* d_V;

      // Rank of process this object lives on.
      int d_rank;

      // Total number of processors.
      int d_size;

      // MPI message tag.
      static const int COMMUNICATE_A;
};

}

#endif
