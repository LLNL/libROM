#ifndef included_IncrementalSVDFastUpdate_h
#define included_IncrementalSVDFastUpdate_h

#include "IncrementalSVD.h"

namespace CAROM {

// A class which embodies the fast update incremental SVD algorithm.
class IncrementalSVDFastUpdate : public IncrementalSVD
{
   public:
      // Constructor.
      IncrementalSVDFastUpdate(
         int dim,
         double redundancy_tol,
         bool skip_redundant,
         int increments_per_time_interval);

      // Destructor.
      ~IncrementalSVDFastUpdate();

      // Increment the SVD with the new state, u_in, at the given time.
      void
      increment(
         const double* u_in,
         double time);

      // Returns the basis vectors for the current time interval, d_U*d_Up, as
      // a Matrix.
      const Matrix*
      getBasis();

   private:
      // Unimplemented default constructor.
      IncrementalSVDFastUpdate();

      // Unimplemented copy constructor.
      IncrementalSVDFastUpdate(
         const IncrementalSVDFastUpdate& other);

      // Unimplemented assignment operator.
      IncrementalSVDFastUpdate&
      operator = (
         const IncrementalSVDFastUpdate& rhs);

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
