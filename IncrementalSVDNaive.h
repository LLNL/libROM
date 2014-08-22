/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The concrete implementation of the incremental SVD algorithm
 *              that is equivalent to but computationally more expensive than
 *              the "fast update" method.
 *
 *****************************************************************************/

#ifndef included_IncrementalSVDNaive_h
#define included_IncrementalSVDNaive_h

#include "IncrementalSVD.h"

namespace CAROM {

// A class which embodies the naive incremental SVD algorithm.
class IncrementalSVDNaive : public IncrementalSVD
{
   public:
      // Constructor.
      IncrementalSVDNaive(
         int dim,
         double redundancy_tol,
         bool skip_redundant,
         int increments_per_time_interval);

      // Destructor.
      ~IncrementalSVDNaive();

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
      IncrementalSVDNaive();

      // Unimplemented copy constructor.
      IncrementalSVDNaive(
         const IncrementalSVDNaive& other);

      // Unimplemented assignment operator.
      IncrementalSVDNaive&
      operator = (
         const IncrementalSVDNaive& rhs);

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
