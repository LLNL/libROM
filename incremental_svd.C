#include "incremental_svd.h"
#include "SLEPcManager.h"
#include <cmath>
#include <assert.h>

extern "C" {
void dgesdd_(char*, int*, int*, double*, int*,
             double*, double*, int*, double*, int*,
             double*, int*, int*, int*);
}

// #define DEBUG

namespace CAROM {

incremental_svd::incremental_svd(
   int* argc,
   char*** argv,
   int dim,
   double epsilon,
   bool skip_redundant) :
   d_dim(dim),
   d_num_increments(0),
   d_epsilon(epsilon),
   d_skip_redundant(skip_redundant),
   d_U(0),
   d_L(0),
   d_S(0),
   d_norm_j(0.0)
{
   // Register a new SLEPc instance, get the rank of this process, and get the
   // number of processors.
   SLEPcManager::getManager()->registerSLEPcInstance(argc, argv);
   MPI_Comm_rank(PETSC_COMM_WORLD, &d_rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &d_size);
}

incremental_svd::~incremental_svd()
{
   // Destroy PETSc data members and finalize SLEPc.
   if (d_U) {
      delete [] d_U;
   }
   if (d_L) {
      delete [] d_L;
   }
   if (d_S) {
      delete [] d_S;
   }
   SLEPcManager::getManager()->unRegisterSLEPcInstance();
}

void
incremental_svd::increment(
   double* u_in)
{
   // If this is the first SVD then build it.  Otherwise add this increment to
   // the system.
   if (d_num_increments == 0) {
      buildInitialSVD(u_in);
   }
   else {
      buildIncrementalSVD(u_in);
   }

#ifdef DEBUG
   if (d_rank == 0) {
      int idx = 0;
      for (int row = 0; row < d_num_increments; ++row) {
         for (int col = 0; col < d_num_increments; ++col) {
            printf("%g  ", d_S[idx++]);
         }
         printf("\n");
      }
      printf("\n");
      idx = 0;
      for (int row = 0; row < d_num_increments; ++row) {
         for (int col = 0; col < d_num_increments; ++col) {
            printf("%g  ", d_L[idx++]);
         }
         printf("\n");
      }
      printf("\n");
   }
   Mat U;
   MatCreate(PETSC_COMM_WORLD, &U);
   MatSetSizes(U, d_dim, PETSC_DECIDE, PETSC_DETERMINE, d_num_increments);
   MatSetType(U, MATDENSE);
   MatSetUp(U);
   int U_row_start;
   MatGetOwnershipRange(U, &U_row_start, PETSC_NULL);
   int* rows = new int [d_dim];
   for (int i = 0; i < d_dim; ++i) {
      rows[i] = i + U_row_start;
   }
   int* cols = new int [d_num_increments];
   for (int i = 0; i < d_num_increments; ++i) {
      cols[i] = i;
   }
   MatSetValues(U, d_dim, rows, d_num_increments, cols, d_U, INSERT_VALUES);
   MatAssemblyBegin(U, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(U, MAT_FINAL_ASSEMBLY);
   MatView(U, PETSC_VIEWER_STDOUT_WORLD);
   MatDestroy(&U);
   if (d_rank == 0) {
      printf("============================================================\n");
   }
#endif
}

void
incremental_svd::buildInitialSVD(
   double* u)
{
   int zero = 0;

   // Build d_S.
   d_S = new double [1];
   double norm_u = 0.0;
   double tmp = 0.0;
   for (int i = 0; i < d_dim; ++i) {
      tmp += u[i]*u[i];
   }
   MPI_Allreduce(&tmp, &norm_u, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
   norm_u = sqrt(norm_u);
   d_S[0] = norm_u;

   // Build d_L.
   d_L = new double [1];
   d_L[0] = 1.0;

   // Build d_U.
   d_U = new double [d_dim];
   for (int i = 0; i < d_dim; ++i) {
      d_U[i] = u[i]/norm_u;
   }

   // We now have another increment.
   ++d_num_increments;
}

void
incremental_svd::buildIncrementalSVD(
   double* u)
{
   // Compute j, P and the norm of j.
   double* j;
   double* P;
   compute_J_P_normJ(u, j, P);

   // l = P' * u
   double* l;
   Pt_Times_u(P, u, l);

   // Compute cm = u.u and see if this increment is new.
   double cm = 0.0;
   double tmp = 0.0;
   for (int i = 0; i < d_dim; ++i) {
      tmp += u[i]*u[i];
   }
   MPI_Allreduce(&tmp, &cm, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

   double eps_squared = d_epsilon*d_epsilon;
   bool is_new_increment = ((cm > eps_squared) &&
      ((d_norm_j*d_norm_j)/(eps_squared+cm)) > eps_squared);

   // If this increment is not new and we are skipping redundant increments
   // then clean up and return.
   double norm_j_for_constructQ = d_norm_j;
   if (!is_new_increment) {
      if (d_skip_redundant) {
         delete [] j;
         delete [] l;
         return;
      }
      else {
         norm_j_for_constructQ = 0.0;
      }
   }

   // Create Q.
   double* Q;
   constructQ(Q, l, norm_j_for_constructQ);

   // Done with l.
   delete [] l;

   // Now get the singular value decomposition of Q.
   double* A;
   double* sigma;
   double* B;
   svd(Q, A, sigma, B);

   // Done with Q.
   delete [] Q;

   // At this point we either have a new or a redundant increment and we are
   // not skipping redundant increments.
   if (!is_new_increment && !d_skip_redundant) {
      // This increment is not new and we are not skipping redundant
      // increments.
      addRedundantIncrement(A, sigma);
      delete [] sigma;
   }
   else if (is_new_increment) {
      // This increment is new.
      addNewIncrement(j, A, sigma);
   }

   // Clean up.
   delete [] P;
   delete [] j;
   delete [] A;
   delete [] B;
}

void
incremental_svd::compute_J_P_normJ(
   double* u,
   double*& j,
   double*& P)
{
   // j = u
   j = new double [d_dim];
   memcpy(j, u, d_dim*sizeof(double));

   // P = d_U * d_L
   d_U_Times_d_L(P);

   // Use modified Gram-Schmidt orthogonalization to modify j.
   orthogonalizeJAndComputeNorm(j, P);
}

void
incremental_svd::orthogonalizeJAndComputeNorm(
   double* j,
   double* P)
{
   double tmp;
   for (int col = 0; col < d_num_increments; ++col) {
      double factor = 0.0;
      tmp = 0.0;
      int Pidx = col;
      for (int i = 0; i < d_dim; ++i) {
         tmp += P[Pidx]*j[i];
         Pidx += d_num_increments;
      }
      MPI_Allreduce(&tmp, &factor, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
      
      Pidx = col;
      for (int i = 0; i < d_dim; ++i) {
         j[i] -= factor*P[Pidx];
         Pidx += d_num_increments;
      }
   }

   // Normalize j.
   tmp = 0.0;
   for (int i = 0; i < d_dim; ++i) {
      tmp += j[i]*j[i];
   }
   MPI_Allreduce(&tmp, &d_norm_j, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
   d_norm_j = sqrt(d_norm_j);
   for (int i = 0; i < d_dim; ++i) {
      j[i] /= d_norm_j;
   }
}

void
incremental_svd::constructQ(
   double*& Q,
   double* l,
   double norm_j)
{
   // Create Q.
   Q = new double [(d_num_increments+1)*(d_num_increments+1)];

   // Fill Q in column major order.
   int q_idx = 0;
   int s_idx = 0;
   for (int row = 0; row < d_num_increments; ++row) {
      q_idx = row;
      for (int col = 0; col < d_num_increments; ++col) {
         Q[q_idx] = d_S[s_idx++];
         q_idx += d_num_increments+1;
      }
      Q[q_idx] = l[row];
   }
   q_idx = d_num_increments;
   for (int col = 0; col < d_num_increments; ++col) {
      Q[q_idx] = 0.0;
      q_idx += d_num_increments+1;
   }
   Q[q_idx] = norm_j;
}

void
incremental_svd::svd(
   double* A,
   double*& U,
   double*& S,
   double*& V)
{
   // Construct U, S, and V.
   int mat_size = (d_num_increments+1)*(d_num_increments+1);
   U = new double [mat_size];
   S = new double [mat_size];
   V = new double [mat_size];
   for (int i = 0; i < mat_size; ++i) {
      S[i] = 0.0;
   }

   // Use lapack's dgesdd_ Fortran function to perform the svd.  As this is
   // Fortran A and all the computed matrices are in column major order.
   double* sigma = new double [d_num_increments+1];
   char jobz = 'A';
   int m = d_num_increments+1;
   int n = d_num_increments+1;
   int lda = d_num_increments+1;
   int ldu = d_num_increments+1;
   int ldv = d_num_increments+1;
   int lwork = m*(4*m + 7);
   double* work = new double [lwork];
   int iwork[8*m];
   int info;
   dgesdd_(&jobz, &m, &n, A, &lda,
           sigma, U, &ldu, V, &ldv,
           work, &lwork, iwork, &info);
   delete [] work;

   // Place sigma into S.
   for (int i = 0; i < d_num_increments+1; ++i) {
      S[i*(d_num_increments+1)+i] = sigma[i];
   }
   delete [] sigma;

   // U and V are in column major order.  Convert them to row major order.
   for (int row = 0; row < d_num_increments+1; ++row) {
      for (int col = row+1; col < d_num_increments+1; ++col) {
         double tmp = U[row*(d_num_increments+1)+col];
         U[row*(d_num_increments+1)+col] = U[col*(d_num_increments+1)+row];
         U[col*(d_num_increments+1)+row] = tmp;
         tmp = V[row*(d_num_increments+1)+col];
         V[row*(d_num_increments+1)+col] = V[col*(d_num_increments+1)+row];
         V[col*(d_num_increments+1)+row] = tmp;
      }
   }
}

void
incremental_svd::addRedundantIncrement(
   double* A,
   double* sigma)
{
   // Chop a row and a column off of A to form Amod.  Also form d_S by chopping
   // a row and a column off of sigma.
   double* Amod = new double [d_num_increments*d_num_increments];
   int lhs_idx = 0;
   int rhs_idx = 0;
   for (int row = 0; row < d_num_increments; ++row){
      for (int col = 0; col < d_num_increments; ++col) {
         Amod[lhs_idx] = A[rhs_idx];
         d_S[lhs_idx] = sigma[rhs_idx];
         ++lhs_idx;
         ++rhs_idx;
      }
      ++rhs_idx;
   }

   // Multiply d_L and Amod and put result into d_L.
   double* L_times_Amod;
   MatTimesMat(d_L, Amod, d_num_increments, L_times_Amod);
   delete [] d_L;
   d_L = L_times_Amod;

   // Clean up.
   delete [] Amod;
}

void
incremental_svd::addNewIncrement(
   double* j,
   double* A,
   double* sigma)
{
   // Add j as a new column of d_U.
   double* newU = new double [d_dim*(d_num_increments+1)];
   int lhs_idx = 0;
   int rhs_idx = 0;
   for (int row = 0; row < d_dim; ++row) {
      for (int col = 0; col < d_num_increments; ++col) {
         newU[lhs_idx++] = d_U[rhs_idx++];
      }
      newU[lhs_idx++] = j[row];
   }
   delete [] d_U;
   d_U = newU;

   // Add another row and column to d_L.  Only the last value in the new
   // row/column is non-zero and it is 1.
   double* newL = new double [(d_num_increments+1)*(d_num_increments+1)];
   lhs_idx = 0;
   rhs_idx = 0;
   for (int row = 0; row < d_num_increments; ++row) {
      for (int col = 0; col < d_num_increments; ++col) {
         newL[lhs_idx++] = d_L[rhs_idx++];
      }
      newL[lhs_idx++] = 0.0;
   }
   for (int col = 0; col < d_num_increments; ++col) {
      newL[lhs_idx++] = 0.0;
   }
   newL[lhs_idx] = 1.0;
   delete [] d_L;
   MatTimesMat(newL, A, d_num_increments+1, d_L);
   delete [] newL;

   // d_S = sigma.
   delete [] d_S;
   d_S = sigma;

   // We now have another increment.
   ++d_num_increments;
}

void
incremental_svd::d_U_Times_d_L(
   double*& P) const
{
   // Create the result.
   P = new double [d_dim * d_num_increments];

   // Construct the product.
   for (int L_col = 0; L_col < d_num_increments; ++L_col) {
      int U_idx = 0;
      int P_idx = L_col;
      for (int U_row = 0; U_row < d_dim; ++U_row) {
         int L_idx = L_col;
         P[P_idx] = 0.0;
         for (int entry = 0; entry < d_num_increments; ++entry) {
            P[P_idx] += d_U[U_idx]*d_L[L_idx];
            ++U_idx;
            L_idx += d_num_increments;
         }
         P_idx += d_num_increments;
      }
   }
}

void
incremental_svd::Pt_Times_u(
   double* P,
   double* u,
   double*& l)
{
   // Create this processor's contribution to the result.
   double* tmp = new double [d_num_increments];
   int l_idx = 0;
   for (int P_col = 0; P_col < d_num_increments; ++P_col) {
      int P_idx = P_col;
      int u_idx = 0;
      tmp[l_idx] = 0.0;
      for (int entry = 0; entry < d_dim; ++entry) {
         tmp[l_idx] += P[P_idx]*u[u_idx];
         P_idx += d_num_increments;
         ++u_idx;
      }
      ++l_idx;
   }

   // Create the final result.
   l = new double [d_num_increments];

   // Sum all processors' contributions into the final result.
   MPI_Allreduce(tmp, l, d_num_increments, MPI_DOUBLE,
                 MPI_SUM, PETSC_COMM_WORLD);

   // Clean up.
   delete [] tmp;
}

void
incremental_svd::MatTimesMat(
   const double* A,
   const double* B,
   int size,
   double*& result)
{
   // Construct result.
   result = new double [size*size];

   for (int B_col = 0; B_col < size; ++B_col) {
      int A_idx = 0;
      int result_idx = B_col;
      for (int A_row = 0; A_row < size; ++A_row) {
         int B_idx = B_col;
         result[result_idx] = 0.0;
         for (int entry = 0; entry < size; ++entry) {
            result[result_idx] += A[A_idx]*B[B_idx];
            ++A_idx;
            B_idx += size;
         }
         result_idx += size;
      }
   }
}

}
