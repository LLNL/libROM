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

/**
 * A class which embodies the naive incremental SVD algorithm.
 */
class IncrementalSVDNaive : public IncrementalSVD
{
   public:
      /**
       * @brief Constructor.
       *
       * @pre dim > 0
       * @pre redundancy_tol > 0.0
       * @pre samples_per_time_interval > 0
       *
       * @param[in] dim The dimension of the system on this processor.
       * @param[in] redundancy_tol Tolerance to determine if a sample is
       *                           redundant or not.
       * @param[in] skip_redundant If true skip redundant samples.
       * @param[in] samples_per_time_interval The number of samples to be
       *                                      collected for each time interval.
       * @param[in] debug_rom If true results of algorithm will be printed to
       *                      facilitate debugging.
       */
      IncrementalSVDNaive(
         int dim,
         double redundancy_tol,
         bool skip_redundant,
         int samples_per_time_interval,
         bool debug_rom = false);

      /**
       * @brief Destructor.
       */
      ~IncrementalSVDNaive();

      /**
       * @brief Returns the basis vectors for the current time interval as a
       * Matrix.
       *
       * @return The basis vectors for the current time interval.
       */
      const Matrix*
      getBasis();

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      IncrementalSVDNaive();

      /**
       * @brief Unimplemented copy constructor.
       */
      IncrementalSVDNaive(
         const IncrementalSVDNaive& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      IncrementalSVDNaive&
      operator = (
         const IncrementalSVDNaive& rhs);

      /**
       * @brief Constructs the first svd.
       *
       * @pre u != 0
       * @pre time >= 0.0
       *
       * @param[in] u The first state.
       * @param[in] time The simulation time for the first state.
       */
      virtual
      void
      buildInitialSVD(
         const double* u,
         double time);

      /**
       * @brief Adds the new sampled state vector, u, to the system.
       *
       * @pre u != 0
       *
       * @param[in] u The new state.
       */
      virtual
      void
      buildIncrementalSVD(
         const double* u);

      /**
       * Add a redundant sample to the svd.
       *
       * @pre A != 0
       * @pre sigma != 0
       *
       * @param[in] A The left singular vectors.
       * @param[in] sigma The singular values.
       */
      void
      addRedundantSample(
         const Matrix* A,
         const Matrix* sigma);

      /**
       * @brief Add a new, unique sample to the svd.
       *
       * @pre j != 0
       * @pre A != 0
       * @pre sigma != 0
       *
       * @param[in] j The new column of d_U.
       * @param[in] A The left singular vectors.
       * @param[in] sigma The singular values.
       */
      void
      addNewSample(
         const Vector* j,
         const Matrix* A,
         Matrix* sigma);

      /**
       * @brief Computes and returns the orthogonality of d_U.
       *
       * @return The orthogonality of d_U.
       */
      double
      checkOrthogonality();

      /**
       * @brief The matrix U distributed across all processors.
       *
       * Each processor's d_U is the part of the distributed matrix local to
       * that processor.
       */
      Matrix* d_U;
};

}

#endif
