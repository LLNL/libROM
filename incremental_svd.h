#include "petscmat.h"
#include "petscvec.h"

namespace CAROM {

// A class which embodies the incremental SVD algorithm.  The API is
// intentionally small.  One may increment the SVD, get the tolerance used
// to determine if an increment is new, and get the dimension of the system.
class incremental_svd
{
   public:
      // Constructor.
      incremental_svd(
         int* argc,
         char*** argv,
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
      getEpsilon()
      {
         return d_epsilon;
      }

      // Returns the dimension of the system.
      int
      getDim()
      {
         return d_dim;
      }

      // Returns the model parameters, d_U*d_L.
      Mat
      getModel()
      {
         Mat result;
         MatMatMult(d_U, d_L, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result);
         return result;
      }

   private:
      friend class incremental_svd_time_stepper;

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
         Vec u);

      // Increments the svd given the state vector u.
      void
      buildIncrementalSVD(
         Vec u);

      // Compute J, P, and the norm of J given u.
      void
      compute_J_P_normJ(
         Vec u,
         Vec& j,
         Mat& P,
         double& norm_j);

      // Helper function which computes the norm of j given u.
      double
      computeNormJ(
         Vec u);

      // Use modified Gram-Schmidt orthogonalization to modify j then compute
      // its norm.
      void
      orthogonalizeJAndComputeNorm(
         Vec j,
         Mat P,
         double& norm_j);

      // Construct the Q matrix which will be passed to svd.
      void
      constructQ(
         Mat& Q,
         Vec l,
         double norm_j);

      // Given a matrix, A, returns the 3 components of the singular value
      // decomposition.
      void
      svd(
         Mat A,
         Mat& U,
         Mat& S,
         Mat& V);

      // Add a redundant increment to the svd.
      void
      addRedundantIncrement(
         Mat A,
         Mat sigma);

      // Add a new, unique increment to the svd.
      void
      addNewIncrement(
         Vec j,
         Mat A,
         Mat sigma,
         int d_U_row_start,
         int d_U_row_end);

      // Dimension of the system.
      int d_dim;

      // Number of increments stored.
      int d_num_increments;

      // Error tolerance.
      double d_epsilon;

      // If true, skip redundant increments.
      bool d_skip_redundant;

      Mat d_U;
      Mat d_L;
      Mat d_S;

      // Rank of process this object lives on.
      int d_rank;
};

}
