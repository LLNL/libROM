/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2015 Lawrence Livermore National Security, LLC
 * Description: The abstract incremental SVD algorithm defines algorithm
 *              interface.
 *
 *****************************************************************************/

#ifndef included_IncrementalSVD_h
#define included_IncrementalSVD_h

#include "SVD.h"

namespace CAROM {

/**
 * IncrementalSVD is an abstract class defining the internal API of the
 * incremental SVD algorithm.
 */
  class IncrementalSVD : public SVD
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
       * @param[in] debug_algorithm If true results of algorithm will be
       *                            printed to facilitate debugging.
       */
      IncrementalSVD(
         int dim,
         double redundancy_tol,
         bool skip_redundant,
         int samples_per_time_interval,
         bool debug_algorithm = false);

      /**
       * @brief Destructor.
       */
      virtual
      ~IncrementalSVD();

      /**
       * @brief Sample new state, u_in, at the given time.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] time The simulation time for the state.
       *
       * @return True if the sampling was successful.
       */
      virtual
      bool
      takeSample(
         const double* u_in,
         double time);

      /**
       * @brief Returns the basis vectors for the current time interval.
       *
       * @return The basis vectors for the current time interval.
       */
      virtual
      const Matrix*
      getBasis();

   protected:
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
         double time) = 0;

      /**
       * @brief Adds the new sampled the state vector, u, to the system.
       *
       * @pre u != 0
       *
       * @param[in] u The new state.
       *
       * @return True if building the incremental svd was successful.
       */
      virtual
      bool
      buildIncrementalSVD(
         const double* u);

      /**
       * @brief Computes the current basis vectors.
       */
      virtual
      void
      computeBasis() = 0;

      /**
       * @brief Construct the matrix Q whose svd is needed.
       *
       * @pre l != 0
       * @pre l.dim() == numSamples()
       *
       * @param[out] Q The matrix to be constructed. [d_S,l; 0,k]
       * @param[in] l The last column of Q.
       * @param[in] k The lower right element of Q.
       */
      void
      constructQ(
         double*& Q,
         const Vector* l,
         double k);

      /**
       * @brief Given a matrix, A, returns 2 of the 3 components of its
       * singular value decomposition.
       *
       * The right singular vectors are not needed and therefore not returned.
       *
       * @pre A != 0
       *
       * @param[in] A The matrix whose svd is needed.
       * @param[out] U The left singular vectors of A.
       * @param[out] S The singular values of A.
       *
       * @return True if the svd succeeded.
       */
      bool
      svd(
         double* A,
         Matrix*& U,
         Matrix*& S);

      /**
       * Add a redundant sample to the svd.
       *
       * @pre A != 0
       * @pre sigma != 0
       *
       * @param[in] A The left singular vectors.
       * @param[in] sigma The singular values.
       */
      virtual
      void
      addRedundantSample(
         const Matrix* A,
         const Matrix* sigma) = 0;

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
      virtual
      void
      addNewSample(
         const Vector* j,
         const Matrix* A,
         Matrix* sigma) = 0;

      /**
       * @brief The number of samples stored.
       *
       * @return The number of samples stored.
       */
      int
      numSamples()
      {
         return d_num_samples;
      }

      /**
       * @brief Computes and returns the orthogonality of m.
       *
       * @pre m != 0
       *
       * @param[in] m The matrix to check.
       *
       * @return The orthogonality of m.
       */
      double
      checkOrthogonality(
        const Matrix* m);

      /**
       * @brief Reorthogonalizes m.
       *
       * @pre m != 0
       *
       * @param[in,out] m The matrix to reorthogonalize.
       */
      void
      reOrthogonalize(
         Matrix* m);

      /**
       * @brief Tolerance to determine if a sample is redundant or not.
       */
      double d_redundancy_tol;

      /**
       * @brief If true, skip redundant samples.
       */
      bool d_skip_redundant;

      /**
       * @brief The matrix S.
       *
       * S is not distributed and the entire matrix exists on each processor.
       */
      Matrix* d_S;

      /**
       * @brief Total number of processors.
       */
      int d_size;

      /**
       * @brief Dimension of the system on each processor.
       */
      std::vector<int> d_proc_dims;

      /**
       * @brief The total dimension of the system.
       */
      long int d_total_dim;

      /**
       * @brief MPI message tag.
       */
      static const int COMMUNICATE_U;

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      IncrementalSVD();

      /**
       * @brief Unimplemented copy constructor.
       */
      IncrementalSVD(
         const IncrementalSVD& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      IncrementalSVD&
      operator = (
         const IncrementalSVD& rhs);
};

}

#endif
