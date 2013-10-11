#include "matrix_utils.h"

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
         bool skip_redundant);

      // Destructor.
      ~incremental_svd();

      // Increment the SVD with the new state, u_in.
      void
      increment(
         double* u_in);

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

      // Returns the number of processors being run on.
      int
      getSize() const
      {
         return d_size;
      }

      // Returns the model parameters, d_U*d_L as a C array.
      double*
      getModel() const
      {
         return DistributedMatLocalMatMult(d_U, d_dim, d_num_increments,
                                           d_L, d_num_increments,
                                           d_num_increments);
      }

      double
      getNormJ() const
      {
         return d_norm_j;
      }

   private:
      // Unimplemented default constructor.
      incremental_svd();

      // Unimplemented copy constructor.
      incremental_svd(const incremental_svd& other);

      // Unimplemented assignment operator.
      incremental_svd&
      operator = (const incremental_svd& rhs);

      // Constructs the first svd.
      void
      buildInitialSVD(
         double* u);

      // Increments the svd given the state vector u.
      void
      buildIncrementalSVD(
         double* u);

      // Compute J, P, and the norm of J given u.
      void
      compute_J_P_normJ(
         double* u,
         double*& j,
         double*& P);

      // Use modified Gram-Schmidt orthogonalization to modify j then compute
      // its norm.
      void
      orthogonalizeJAndComputeNorm(
         double* j,
         double* P);

      // Construct the Q matrix which will be passed to svd.
      void
      constructQ(
         double*& Q,
         double* l,
         double norm_j);

      // Given a matrix, A, returns the 3 components of the singular value
      // decomposition.
      void
      svd(
         double* A,
         double*& U,
         double*& S,
         double*& V);

      // Add a redundant increment to the svd.
      void
      addRedundantIncrement(
         double* A,
         double* sigma);

      // Add a new, unique increment to the svd.
      void
      addNewIncrement(
         double* j,
         double* A,
         double* sigma);

      // Dimension of the system.
      int d_dim;

      // Number of increments stored.
      int d_num_increments;

      // Error tolerance.
      double d_epsilon;

      // If true, skip redundant increments.
      bool d_skip_redundant;

      // The matrix U distributed across all processors.
      double* d_U;

      // The matrix L.  L is not distributed and the entire matrix exists on
      // each processor.
      double* d_L;

      // The matrix S.  S is not distributed and the entire matrix exists on
      // each processor.
      double* d_S;

      // Rank of process this object lives on.
      int d_rank;

      // Total number of processors.
      int d_size;

      // Value of norm of j cached by increment.
      double d_norm_j;

      // MPI message tag.
      static const int COMMUNICATE_U;
};

}
