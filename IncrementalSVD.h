/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The abstract incremental SVD algorithm defines algorithm
 *              interface.
 *
 *****************************************************************************/

#ifndef included_IncrementalSVD_h
#define included_IncrementalSVD_h

#include "Matrix.h"
#include <vector>

namespace CAROM {

/**
 * IncrementalSVD is an abstract class defining the API of the incremental SVD
 * algorithm.  The API is intentionally small.  One may take a sample, get the
 * dimension of the system, get the basis vectors, get the number of time
 * intervals and their start times, and determine if the current time interval
 * is full.
 */
class IncrementalSVD
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
      IncrementalSVD(
         int dim,
         double redundancy_tol,
         bool skip_redundant,
         int samples_per_time_interval,
         bool debug_rom = false);

      /**
       * @brief Destructor.
       */
      virtual
      ~IncrementalSVD();

      /**
       * @brief Sample the new state, u_in, at the given time.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] time The simulation time for the state.
       */
      virtual
      void
      takeSample(
         const double* u_in,
         double time) = 0;

      /**
       * @brief Returns the dimension of the system on this processor.
       *
       * @return The dimension of the system on this processor.
       */
      int
      getDim() const
      {
         return d_dim;
      }

      /**
       * @brief Returns the basis vectors for the current time interval as a
       * Matrix.
       *
       * @return The basis vectors for the current time interval.
       */
      virtual
      const Matrix*
      getBasis() = 0;

      /**
       * @brief Returns the number of time intervals on which different sets
       * of basis vectors are defined.
       *
       * @return The number of time intervals for which basis vectors exist.
       */
      int
      getNumBasisTimeIntervals() const
      {
         return static_cast<int>(d_time_interval_start_times.size());
      }

      /**
       * @brief Returns the start time for the requested time interval.
       *
       * @pre 0 <= which_interval
       * @pre which_interval < getNumBasisTimeIntervals()
       *
       * @param[in] which_interval Time interval whose start time is needed.
       *
       * @return The start time of the requested time interval.
       */
      double
      getBasisIntervalStartTime(
         int which_interval) const
      {
         CAROM_ASSERT(0 <= which_interval);
         CAROM_ASSERT(which_interval < getNumBasisTimeIntervals());
         return d_time_interval_start_times[which_interval];
      }

      /**
       * @brief Returns true if the next state will result in a new time
       * interval.
       *
       * @return True if all samples have been taken for the current time
       * interval.
       */
      bool
      isNewTimeInterval() const
      {
         return (d_num_samples == 0) ||
                (d_num_samples >= d_samples_per_time_interval);
      }

   protected:
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
       */
      void
      svd(
         double* A,
         Matrix*& U,
         Matrix*& S);

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
       * @brief Dimension of the system on this processor.
       */
      int d_dim;

      /**
       * @brief Number of samples stored.
       */
      int d_num_samples;

      /**
       * @brief Tolerance to determine if a sample is redundant or not.
       */
      double d_redundancy_tol;

      /**
       * @brief If true, skip redundant samples.
       */
      bool d_skip_redundant;

      /**
       * @brief The number of samples to be collected for each time interval.
       */
      int d_samples_per_time_interval;

      /**
       * @brief The matrix S.
       *
       * S is not distributed and the entire matrix exists on each processor.
       */
      Matrix* d_S;

      /**
       * @brief The basis vectors for the current time interval distributed
       * across all processors.
       *
       * d_basis is the part of the distributed basis vector local to the
       * processor owning this object.
       */
      Matrix* d_basis;

      /**
       * @brief The simulation time at which each time interval starts.
       */
      std::vector<double> d_time_interval_start_times;

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
       * @brief Flag to indicate if results of algorithm should be printed for
       * debugging purposes.
       */
      bool d_debug_rom;

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
