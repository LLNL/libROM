#ifndef included_incremental_svd_fast_update_h
#define included_incremental_svd_fast_update_h

#include "incremental_svd.h"

namespace CAROM {

// A class which embodies the fast update incremental SVD algorithm.
class incremental_svd_fast_update : public incremental_svd
{
   public:
      // Constructor.
      incremental_svd_fast_update(
         int dim,
         double epsilon,
         bool skip_redundant,
         int increments_per_time_interval);

      // Destructor.
      ~incremental_svd_fast_update();

      // Increment the SVD with the new state, u_in, at the given time.
      void
      increment(
         const double* u_in,
         double time);

      // Returns the basis vectors for the given time, d_U*d_Up, as a Matrix.
      const Matrix*
      getBasis(
         double time);

   private:
      // Unimplemented default constructor.
      incremental_svd_fast_update();

      // Unimplemented copy constructor.
      incremental_svd_fast_update(
         const incremental_svd_fast_update& other);

      // Unimplemented assignment operator.
      incremental_svd_fast_update&
      operator = (
         const incremental_svd_fast_update& rhs);

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
      checkUOrthogonality();

      // Returns the orthogonality of d_Up.
      double
      checkUpOrthogonality();

      // Reorthogonalizes d_U.
      void
      reOrthogonalizeU();

      // Reorthogonalizes d_Up.
      void
      reOrthogonalizeUp();

      // The matrix U distributed across all processors.  Each processor's d_U
      // is the part of the distributed matrix local to that processor.
      Matrix* d_U;

      // The matrix U'.  U' is not distributed and the entire matrix exists on
      // each processor.
      Matrix* d_Up;
};

}

#endif
