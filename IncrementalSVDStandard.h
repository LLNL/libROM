/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2015 Lawrence Livermore National Security, LLC
 * Description: The concrete implementation of the incremental SVD algorithm
 *              that is equivalent to but computationally more expensive than
 *              the "fast update" method.
 *
 *****************************************************************************/

#ifndef included_IncrementalSVDStandard_h
#define included_IncrementalSVDStandard_h

#include "IncrementalSVD.h"

namespace CAROM {

/**
 * A class which embodies the standard incremental SVD algorithm.
 */
class IncrementalSVDStandard : public IncrementalSVD
{
   public:
      /**
       * @brief Constructor.
       *
       * @pre dim > 0
       * @pre linearity_tol > 0.0
       * @pre samples_per_time_interval > 0
       *
       * @param[in] dim The dimension of the system on this processor.
       * @param[in] linearity_tol Tolerance to determine whether or not a
       *                          sample is linearly dependent.
       * @param[in] skip_linearly_dependent If true skip linearly dependent
       *                                    samples.
       * @param[in] samples_per_time_interval The number of samples to be
       *                                      collected for each time interval.
       * @param[in] debug_algorithm If true results of the algorithm will be
       *                            printed to facilitate debugging.
       */
      IncrementalSVDStandard(
         int dim,
         double linearity_tol,
         bool skip_linearly_dependent,
         int samples_per_time_interval,
         bool debug_algorithm = false);

      /**
       * @brief Destructor.
       */
      ~IncrementalSVDStandard();

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      IncrementalSVDStandard();

      /**
       * @brief Unimplemented copy constructor.
       */
      IncrementalSVDStandard(
         const IncrementalSVDStandard& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      IncrementalSVDStandard&
      operator = (
         const IncrementalSVDStandard& rhs);

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
       * @brief Computes the current basis vectors.
       */
      virtual
      void
      computeBasis();

      /**
       * Add a linearly dependent sample to the svd.
       *
       * @pre A != 0
       * @pre sigma != 0
       *
       * @param[in] A The left singular vectors.
       * @param[in] sigma The singular values.
       */
      void
      addLinearlyDependentSample(
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
       * @brief The matrix U distributed across all processors.
       *
       * Each processor's d_U is the part of the distributed matrix local to
       * that processor.
       */
      Matrix* d_U;
};

}

#endif
