#ifndef included_incremental_svd_h
#define included_incremental_svd_h

#include "HDFDatabase.h"
#include "matrix.h"
#include <vector>

namespace CAROM {

// An abstract class which embodies the incremental SVD algorithm.  The API is
// intentionally small.  One may increment the SVD, get the tolerance used
// to determine if an increment is new, and get the dimension of the system.
class incremental_svd
{
   public:
      // Constructor.
      incremental_svd(
         int dim,
         double epsilon,
         bool skip_redundant,
         int increments_per_time_interval);

      // Destructor.
      virtual ~incremental_svd();

      // Increment the SVD with the new state, u_in, at the given time.
      virtual
      void
      increment(
         const double* u_in,
         double time) = 0;

      // Returns the tolerance used to determine if an increment is new.
      double
      getEpsilon() const
      {
         return d_epsilon;
      }

      // Returns the dimension of the system.
      int
      getDim() const
      {
         return d_dim;
      }

      // Returns the rank of the processor being run on.
      int
      getRank() const
      {
         return d_rank;
      }

      // Returns the number of processors being run on.
      int
      getSize() const
      {
         return d_size;
      }

      // Returns the basis vectors for the given time, d_U, as a Matrix.
      virtual
      const Matrix*
      getBasis(
         double time) = 0;

      // Returns the value of the norm of J cached by increment.
      double
      getNormJ() const
      {
         return d_norm_j;
      }

      // Writes the basis vectors to a file with the given base name.
      void
      writeBasis(
         const std::string& base_file_name);

   protected:
      // Construct the Q matrix which will be passed to svd.
      void
      constructQ(
         double*& Q,
         const Vector* l,
         double k);

      // Given a matrix, A, returns the 3 components of the singular value
      // decomposition.
      void
      svd(
         double* A,
         Matrix*& U,
         Matrix*& S);

   protected:

      // Dimension of the system.
      int d_dim;

      // Number of increments stored.
      int d_num_increments;

      // Error tolerance.
      double d_epsilon;

      // If true, skip redundant increments.
      bool d_skip_redundant;

      // The number of increments to be collected for each time interval.
      int d_increments_per_time_interval;

      // The matrix S.  S is not distributed and the entire matrix exists on
      // each processor.
      Matrix* d_S;

      // For each time interval, the basis vectors distributed across all
      // processors.  Each d_basis is the part of the distributed basis vectors
      // local to the processor owning this object.
      std::vector<Matrix*> d_basis;

      // The number of time intervals gathered.
      int d_num_time_intervals;

      // The simulation time at which each time interval starts.
      std::vector<double> d_time_interval_start_times;

      // Rank of process this object lives on.
      int d_rank;

      // Total number of processors.
      int d_size;

      // Value of norm of j cached by increment.
      double d_norm_j;

      // Database to read/write basis vectors from/to.
      HDFDatabase d_database;

      // MPI message tag.
      static const int COMMUNICATE_U;

   private:
      // Unimplemented default constructor.
      incremental_svd();

      // Unimplemented copy constructor.
      incremental_svd(
         const incremental_svd& other);

      // Unimplemented assignment operator.
      incremental_svd&
      operator = (
         const incremental_svd& rhs);
};

}

#endif
