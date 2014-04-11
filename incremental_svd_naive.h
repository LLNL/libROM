#ifndef included_incremental_svd_naive_h
#define included_incremental_svd_naive_h

#include "incremental_svd.h"

namespace CAROM {

// A class which embodies the naive incremental SVD algorithm.
class incremental_svd_naive : public incremental_svd
{
   public:
      // Constructor.
      incremental_svd_naive(
         int dim,
         double epsilon,
         bool skip_redundant,
         int increments_per_time_interval);

      // Destructor.
      ~incremental_svd_naive();

      // Increment the SVD with the new state, u_in, at the given time.
      void
      increment(
         const double* u_in,
         double time);

      // Returns the basis vectors for the current time interval, d_U, as a
      // Matrix.
      const Matrix*
      getBasis();

   private:
      // Unimplemented default constructor.
      incremental_svd_naive();

      // Unimplemented copy constructor.
      incremental_svd_naive(
         const incremental_svd_naive& other);

      // Unimplemented assignment operator.
      incremental_svd_naive&
      operator = (
         const incremental_svd_naive& rhs);

      // Constructs the first svd.
      void
      buildInitialSVD(
         const double* u,
         double time);

      // Increments the svd given the state vector u.
      void
      buildIncrementalSVD(
         const double* u);

      // Add a redundant increment to the svd.
      void
      addRedundantIncrement(
         const Matrix* A,
         const Matrix* sigma);

      // Add a new, unique increment to the svd.
      void
      addNewIncrement(
         const Vector* j,
         const Matrix* A,
         Matrix* sigma);

      // Returns the orthogonality of d_U.
      double
      checkOrthogonality();

      // Reorthogonalizes d_U.
      void
      reOrthogonalize();

      // The matrix U distributed across all processors.  Each processor's d_U
      // is the part of the distributed matrix local to that processor.
      Matrix* d_U;
};

}

#endif
