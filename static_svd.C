#include "static_svd.h"
#include "SLEPcManager.h"
#include "slepc.h"
#include <assert.h>

// #define DEBUG

namespace CAROM {

static_svd::static_svd(
   int* argc,
   char*** argv,
   int dim) :
   d_dim(dim),
   d_state(0),
   d_U(0),
   d_S(0),
   d_V(0)
{
   // Register a new SLEPc instance and get the rank of this process.
   SLEPcManager::getManager()->registerSLEPcInstance(argc, argv);
   MPI_Comm_rank(PETSC_COMM_WORLD, &d_rank);
}

static_svd::~static_svd()
{
   // Destroy PETSc data members and unregister a SLEPc instance.
   if (d_U) {
      MatDestroy(&d_U);
   }
   if (d_S) {
      MatDestroy(&d_S);
   }
   if (d_V) {
      MatDestroy(&d_V);
   }
   SLEPcManager::getManager()->unRegisterSLEPcInstance();
}

void
static_svd::computeSVD()
{
   // Construct A.
   int num_cols = static_cast<int>(d_state.size());
   int* rows = new int[d_dim];
   int start_row = d_dim*d_rank;
   for (int i = 0; i < d_dim; ++i) {
      rows[i] = start_row + i;
   }
   Mat A;
   MatCreate(PETSC_COMM_WORLD, &A);
   MatSetSizes(A, d_dim, PETSC_DECIDE, PETSC_DETERMINE, num_cols);
   MatSetType(A, MATDENSE);
   MatSetUp(A);
   for (int i = 0; i < num_cols; ++i) {
      MatSetValues(A, d_dim, rows, 1, &i, d_state[i], INSERT_VALUES);
   }
   MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
   delete [] rows;
   svd(A);
#ifdef DEBUG
   MatView(d_U, PETSC_VIEWER_STDOUT_WORLD);
   MatView(d_S, PETSC_VIEWER_STDOUT_WORLD);
#endif
}

void
static_svd::svd(
   Mat A)
{
   int num_states = static_cast<int>(d_state.size());

   Vec u, v;
   MatGetVecs(A, &v, &u);

   // Construct d_U.
   MatCreate(PETSC_COMM_WORLD, &d_U);
   MatSetSizes(d_U, d_dim, PETSC_DECIDE, PETSC_DETERMINE, num_states);
   MatSetType(d_U, MATDENSE);
   MatSetUp(d_U);

   // Construct d_S.
   MatCreate(PETSC_COMM_WORLD, &d_S);
   MatSetSizes(d_S, PETSC_DECIDE, PETSC_DECIDE, num_states, num_states);
   MatSetType(d_S, MATDENSE);
   MatSetUp(d_S);

   // Construct V.
   MatCreate(PETSC_COMM_WORLD, &d_V);
   MatSetSizes(d_V, PETSC_DECIDE, PETSC_DECIDE, num_states, num_states);
   MatSetType(d_V, MATDENSE);
   MatSetUp(d_V);

   // Solve for the singular value decomposition.
   SVD svd;
   SVDCreate(PETSC_COMM_WORLD, &svd);
   SVDSetOperator(svd, A);
   SVDSetType(svd, SVDLANCZOS);
   SVDSolve(svd);

   int num_sing_vals;
   SVDGetConverged(svd, &num_sing_vals);
   assert(num_sing_vals == static_cast<int>(d_state.size()));

   // The SLEPc SVD solver does not explicitly export U, S, and V so we
   // need to query the solver for their values and construct them ourself.
   // Note that the solver holds the transpose of V so we need be careful
   // when constructing it.
   int u_row_start, u_row_end;
   VecGetOwnershipRange(u, &u_row_start, &u_row_end);
   int num_u_rows = u_row_end - u_row_start;
   int* u_rows = 0;
   double* uvec = 0;
   if (num_u_rows > 0) {
      u_rows = new int[num_u_rows];
      for (int i = 0; i < num_u_rows; ++i) {
         u_rows[i] = u_row_start + i;
      }
      uvec = new double [num_u_rows];
   }

   int v_row_start, v_row_end;
   VecGetOwnershipRange(v, &v_row_start, &v_row_end);
   int num_v_rows = v_row_end - v_row_start;
   int* v_rows = 0;
   double* vvec = 0;
   if (num_v_rows > 0) {
      v_rows = new int [num_v_rows];
      for (int i = 0; i < num_v_rows; ++i) {
         v_rows[i] = v_row_start + i;
      }
      vvec = new double [num_v_rows];
   }
   int* sigma_rows = new int[num_states];
   double* sigma = new double[num_states];
   for (int i = 0; i < num_states; ++i) {
      sigma[i] = 0.0;
      sigma_rows[i] = i;
   }
   for (int col = 0; col < num_states; ++col) {
      if (num_u_rows > 0) {
         if (num_v_rows > 0) {
            SVDGetSingularTriplet(svd, col, &sigma[col], u, v);
            VecGetValues(v, num_v_rows, v_rows, vvec);
            MatSetValues(d_V, num_v_rows, v_rows, 1, &col, vvec, INSERT_VALUES);
            MatSetValues(d_S, num_states, sigma_rows, 1, &col,
                         sigma, INSERT_VALUES);
         }
         else {
            SVDGetSingularTriplet(svd, col, &sigma[col], u, PETSC_NULL);
         }
         VecGetValues(u, num_u_rows, u_rows, uvec);
         MatSetValues(d_U, num_u_rows, u_rows, 1, &col, uvec, INSERT_VALUES);
      }
      else if (num_v_rows > 0) {
         SVDGetSingularTriplet(svd, col, &sigma[col], PETSC_NULL, v);
         VecGetValues(v, num_v_rows, v_rows, vvec);
         MatSetValues(d_V, num_v_rows, v_rows, 1, &col, vvec, INSERT_VALUES);
         MatSetValues(d_S, num_states, sigma_rows, 1, &col,
                      sigma, INSERT_VALUES);
      }
      sigma[col] = 0.0;
   }
   if (u_rows) {
      delete [] u_rows;
      delete [] uvec;
   }
   if (v_rows) {
      delete [] v_rows;
      delete [] vvec;
   }
   delete [] sigma_rows;
   delete [] sigma;
   SVDDestroy(&svd);
   VecDestroy(&u);
   VecDestroy(&v);

   // Assemble d_U, d_S, and d_V.
   MatAssemblyBegin(d_U, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(d_U, MAT_FINAL_ASSEMBLY);
   MatAssemblyBegin(d_S, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(d_S, MAT_FINAL_ASSEMBLY);
   MatAssemblyBegin(d_V, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(d_V, MAT_FINAL_ASSEMBLY);

   // We want to return the transpose of d_V, not d_V.
   MatTranspose(d_V, MAT_REUSE_MATRIX, &d_V);
}

}
