#ifndef included_StaticSVD_h
#define included_StaticSVD_h

#include "Matrix.h"
#include <vector>

namespace CAROM {

// A class which embodies the static SVD algorithm.  The API is intentionally
// small.  One may collect the states, compute the SVD, and get the dimension
// of the system.
// This implementation is not scalable and is intended primarily as a sanity
// check of the incremental svd algorithm.
class StaticSVD
{
   public:
      // Constructor.
      StaticSVD(
         int dim,
         int increments_per_time_interval);

      // Destructor.
      ~StaticSVD();

      // Collect the new state, u_in.
      void
      collectState(
         double* u_in,
         double time);

      // Returns the basis vectors for the current time interval, d_U.
      const Matrix*
      getBasis();

      // Returns the number of time intervals on which different sets of basis
      // vectors are defined.
      int
      getNumBasisTimeIntervals() const
      {
         return d_num_time_intervals;
      }

      // Returns the start time for the requested time interval.
      double
      getBasisIntervalStartTime(
         int which_interval) const
      {
         CAROM_ASSERT(0 <= which_interval);
         CAROM_ASSERT(which_interval < d_num_time_intervals);
         return d_time_interval_start_times[which_interval];
      }

      // Returns true if the next state will result in a new time interval.
      bool
      isNewTimeInterval() const
      {
         return (d_num_increments == 0) ||
                (d_num_increments >= d_increments_per_time_interval);
      }

   private:
      // Unimplemented default constructor.
      StaticSVD();

      // Unimplemented copy constructor.
      StaticSVD(
         const StaticSVD& other);

      // Unimplemented assignment operator.
      StaticSVD&
      operator = (
         const StaticSVD& rhs);

      // Compute the SVD.
      void
      computeSVD();

      // Given a matrix, A, returns the 3 components of the singular value
      // decomposition.
      void
      svd(
         double* A);

      // Dimension of the system.
      int d_dim;

      // Number of increments stored.
      int d_num_increments;

      // The number of increments to be collected for each time interval.
      int d_increments_per_time_interval;

      // Current state of the system.
      std::vector<double*> d_state;

      // The globalized matrix U.  U is large and each process owns all of U.
      Matrix* d_U;

      // The globalized matrix S.  S is small and each process owns all of S.
      Matrix* d_S;

      // The globalized matrix L.  L is small and each process owns all of L.
      Matrix* d_V;

      // The globalized basis vectors for the current time interval.  The basis
      // vectors are large and each process owns all of the basis vectors.
      Matrix* d_basis;

      // The number of time intervals gathered.
      int d_num_time_intervals;

      // The simulation time at which each time interval starts.
      std::vector<double> d_time_interval_start_times;

      // Rank of process this object lives on.
      int d_rank;

      // Total number of processors.
      int d_size;

      // Flag to indicate if the basis vectors for the current time interval
      // are up to date.
      bool d_this_interval_basis_current;

      // MPI message tag.
      static const int COMMUNICATE_A;
};

}

#endif
