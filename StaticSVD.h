/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: A class implementing the static SVD algorithm.
 *
 *****************************************************************************/

#ifndef included_StaticSVD_h
#define included_StaticSVD_h

#include "Matrix.h"
#include <vector>

namespace CAROM {

/**
 * StaticSVD embodies the static SVD algorithm.  The API is intentionally
 * small.  One may collect the states, compute the SVD, and get the dimension
 * of the system.
 * This algorithm is not scalable and is intended primarily as a sanity check
 * of the incremental svd algorithm.
 */
class StaticSVD
{
   public:
      /**
       * @brief Constructor.
       *
       * @pre dim > 0
       * @pre states_per_time_interval > 0
       *
       * @param[in] dim The dimension of the system distributed to this
       *                processor.
       * @param[in] states_per_time_interval The maximium number of states
       *                                     collected in a time interval.
       */
      StaticSVD(
         int dim,
         int states_per_time_interval);

      /**
       * Destructor.
       */
      ~StaticSVD();

      /**
       * @brief Collect the new state, u_in at supplied time.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The new state.
       * @param[in] time The simulation time of the new state.
       */
      void
      collectState(
         double* u_in,
         double time);

      /**
       * @brief Returns the basis vectors for the current time interval.
       *
       * @pre !d_this_interval_basis_current || d_basis != 0
       *
       * @post d_this_interval_basis_current
       *
       * @return The basis vectors for the current time interval.
       */
      const Matrix*
      getBasis();

      /**
       * @brief Returns the number of time intervals on which different sets
       * of basis vectors are defined.
       *
       * @return The number of time intervals on which there are basis vectors.
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
       * @pre which_interval < static_cast<int>(d_time_interval_start_times.size())
       *
       * @param[in] which_interval The time interval of interest.
       *
       * @return The start time for the requested time interval.
       */
      double
      getBasisIntervalStartTime(
         int which_interval) const
      {
         CAROM_ASSERT(0 <= which_interval);
         CAROM_ASSERT(which_interval < static_cast<int>(d_time_interval_start_times.size()));
         return d_time_interval_start_times[which_interval];
      }

      /**
       * @brief Returns true if the next state will result in a new time
       * interval.
       *
       * @return True if the next state results in the creation of a new time
       *         interval.
       */
      bool
      isNewTimeInterval() const
      {
         return (d_num_states == 0) ||
                (d_num_states >= d_states_per_time_interval);
      }

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      StaticSVD();

      /**
       * @brief Unimplemented copy constructor.
       */
      StaticSVD(
         const StaticSVD& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      StaticSVD&
      operator = (
         const StaticSVD& rhs);

      /**
       * @brief Gathers states from all other processors to form complete
       * state of system and computes the SVD.
       */
      void
      computeSVD();

      /**
       * @brief Preforms the actual SVD via lapack.
       *
       * Given a matrix, A, computes the 3 components of the singular value
       * decomposition.
       *
       * @pre A != 0
       * @pre total_dim > 0
       *
       * @param[in] A The globalized system whose SVD will be computed.
       * @param[in] total_dim The total dimension of the system that has been
       *                      distributed of multiple processors.
       */
      void
      svd(
         double* A,
         int total_dim);

      /**
       * @brief Dimension of the system.
       */
      int d_dim;

      /**
       * @brief Number of states stored for the current time interval.
       */
      int d_num_states;

      /**
       * @brief The maximum number of states to be collected for a time
       * interval.
       */
      int d_states_per_time_interval;

      /**
       * @brief Current state of the system.
       */
      std::vector<double*> d_state;

      /**
       * @brief The globalized matrix U.
       *
       * U is large and each process owns all of U.
       */
      Matrix* d_U;

      /**
       * @brief The globalized matrix S.
       *
       * S is small and each process owns all of S.
       */
      Matrix* d_S;

      /**
       * @brief The globalized matrix L.
       *
       * L is small and each process owns all of L.
       */
      Matrix* d_V;

      /**
       * @brief The globalized basis vectors for the current time interval.
       *
       * The basis vectors are large and each process owns all of the basis
       * vectors.
       */
      Matrix* d_basis;

      /**
       * @brief The simulation time at which each time interval starts.
       */
      std::vector<double> d_time_interval_start_times;

      /**
       * @brief Flag to indicate if the basis vectors for the current time
       * interval are up to date.
       */
      bool d_this_interval_basis_current;

      /**
       * @brief MPI message tag.
       */
      static const int COMMUNICATE_A;
};

}

#endif
