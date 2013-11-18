#ifndef included_incremental_svd_h
#define included_incremental_svd_h

#include "HDFDatabase.h"
#include "matrix.h"
#include "vector.h"
#include <vector>

namespace CAROM {

// A class which embodies the incremental SVD algorithm.  The API is
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
      ~incremental_svd();

      // Increment the SVD with the new state, u_in, at the given time.
      void
      increment(
         double* u_in,
         double time);

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

      // Returns the model parameters for the given time, d_U*d_L, as a Matrix.
      const Matrix*
      getModel(
         double time);

      // Returns the value of the norm of J cached by increment.
      double
      getNormJ() const
      {
         return d_norm_j;
      }

      // Writes the model to a file with the given base name.
      void
      writeModel(
         const std::string& base_file_name);

      // Reads the model from a file with the given base name.
      void
      readModel(
         const std::string& base_file_name);

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

      // Constructs the first svd.
      void
      buildInitialSVD(
         double* u,
         double time);

      // Increments the svd given the state vector u.
      void
      buildIncrementalSVD(
         double* u);

      // Compute J, P, and the norm of J given u.
      void
      compute_J_P_normJ(
         const double* u,
         Vector*& j,
         Matrix*& P);

      // Use modified Gram-Schmidt orthogonalization to modify j then compute
      // its norm.
      void
      orthogonalizeJAndComputeNorm(
         Vector* j,
         const Matrix* P);

      // Construct the Q matrix which will be passed to svd.
      void
      constructQ(
         double*& Q,
         const Vector* l,
         double norm_j);

      // Given a matrix, A, returns the 3 components of the singular value
      // decomposition.
      void
      svd(
         double* A,
         Matrix*& U,
         Matrix*& S,
         Matrix*& V);

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

      // For each time interval, the matrix U distributed across all
      // processors.  Each d_U is the part of the distributed matrix local to
      // the processor owning this object.
      std::vector<Matrix*> d_U;

      // For each time interval, the matrix L.  L is not distributed and the
      // entire matrix exists on each processor.
      std::vector<Matrix*> d_L;

      // For each time interval, the matrix S.  S is not distributed and the
      // entire matrix exists on each processor.
      std::vector<Matrix*> d_S;

      // For each time interval, the model parameters distributed across all
      // processors.  Each d_model is the part of the distributed model
      // parameters local to the processor owning this object.
      std::vector<Matrix*> d_model;

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

      // Database to read/write model parameters from/to.
      HDFDatabase database;

      // MPI message tag.
      static const int COMMUNICATE_U;
};

}

#endif
