#ifndef included_static_svd_h
#define included_static_svd_h

#include "HDFDatabase.h"
#include "matrix.h"
#include <vector>
#include <string.h>

namespace CAROM {

// A class which embodies the static SVD algorithm.  The API is intentionally
// small.  One may collect the states, compute the SVD, and get the dimension
// of the system.
// This implementation is not scalable and is intended primarily as a sanity
// check of the incremental svd algorithm.
class static_svd
{
   public:
      // Constructor.
      static_svd(
         int dim,
         int increments_per_time_interval);

      // Destructor.
      ~static_svd();

      // Collect the new state, u_in.
      void
      collectState(
         double* u_in,
         double time);

      // Returns the dimension of the system.
      int
      getDim() const
      {
         return d_dim;
      }

      // Returns the basis vectors, d_U.
      const Matrix*
      getBasis(
         double time);

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

   private:
      // Unimplemented default constructor.
      static_svd();

      // Unimplemented copy constructor.
      static_svd(
         const static_svd& other);

      // Unimplemented assignment operator.
      static_svd&
      operator = (
         const static_svd& rhs);

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

      // For each time interval, the globalized basis vectors.  The basis is
      // large and each process owns all of the basis.
      std::vector<Matrix*> d_basis;

      // The number of time intervals gathered.
      int d_num_time_intervals;

      // The simulation time at which each time interval starts.
      std::vector<double> d_time_interval_start_times;

      // Rank of process this object lives on.
      int d_rank;

      // Total number of processors.
      int d_size;

      // Database to read/write basis vectors from/to.
      HDFDatabase d_database;

      // MPI message tag.
      static const int COMMUNICATE_A;
};

}

#endif
