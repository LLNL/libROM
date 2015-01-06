/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The concrete implementation of the incremental SVD algorithm
 *              using Matthew Brand's "fast update" method.
 *
 *****************************************************************************/

#ifndef included_IncrementalSVDFastUpdate_h
#define included_IncrementalSVDFastUpdate_h

#include "IncrementalSVD.h"

namespace CAROM {

/**
 * IncrementalSVDFastUpdate implements Brand's fast update incremental SVD
 * algorithm by implementing the pure virtual methods of the IncrementalSVD
 * base class.
 */
class IncrementalSVDFastUpdate : public IncrementalSVD
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
      IncrementalSVDFastUpdate(
         int dim,
         double redundancy_tol,
         bool skip_redundant,
         int samples_per_time_interval,
         bool debug_rom = false);

      /**
       * @brief Destructor.
       */
      ~IncrementalSVDFastUpdate();

      /**
       * @brief Sample new state, u_in, at the given time.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] time The simulation time for the state.
       */
      void
      takeSample(
         const double* u_in,
         double time);

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
      IncrementalSVDFastUpdate();

      /**
       * @brief Unimplemented copy constructor.
       */
      IncrementalSVDFastUpdate(
         const IncrementalSVDFastUpdate& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      IncrementalSVDFastUpdate&
      operator = (
         const IncrementalSVDFastUpdate& rhs);

      /**
       * @brief Constructs the first svd.
       *
       * @pre u != 0
       * @pre time >= 0.0
       *
       * @param[in] u The first state.
       * @param[in] time The simulation time for the first state.
       */
      void
      buildInitialSVD(
         const double* u,
         double time);

      /**
       * @brief Adds the new sampled the state vector, u, to the system.
       *
       * @pre u != 0
       *
       * @param[in] u The new state.
       */
      void
      buildIncrementalSVD(
         const double* u);

      /**
       * @brief Add a redundant sample to the svd.
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
      checkUOrthogonality();

      /**
       * @brief Computes and returns the orthogonality of d_Up.
       *
       * @return The orthogonality of d_U.
       */
      double
      checkUpOrthogonality();

      /**
       * @brief Reorthogonalizes m.
       *
       * @pre m != 0
       *
       * @param[in/out] The matrix to reorthogonalize.
       */
      void
      reOrthogonalize(
         Matrix* m);

      /**
       * @brief The matrix U distributed across all processors.
       *
       * Each processor's d_U is the part of the distributed matrix local to
       * that processor.
       */
      Matrix* d_U;

      /**
       * @brief The matrix U'.
       *
       * U' is not distributed and the entire matrix exists on each processor.
       */
      Matrix* d_Up;
};

}

#endif
