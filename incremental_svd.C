#include "incremental_svd.h"
#include "SLEPcManager.h"
#include "slepc.h"
#include <cmath>
#include <assert.h>

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
   d_S(0)
{
   // Register a new SLEPc instance and get the rank of this process.
   SLEPcManager::getManager()->registerSLEPcInstance(argc, argv);
   MPI_Comm_rank(PETSC_COMM_WORLD, &d_rank);
}

incremental_svd::~incremental_svd()
{
   // Destroy PETSc data members and finalize SLEPc.
   if (d_U) {
      MatDestroy(&d_U);
   }
   if (d_L) {
      MatDestroy(&d_L);
   }
   if (d_S) {
      MatDestroy(&d_S);
   }
   SLEPcManager::getManager()->unRegisterSLEPcInstance();
}

void
incremental_svd::increment(
   double* u_in)
{
   // Convert representation free C array of doubles into a PETSc vector.
   Vec u;
   VecCreate(PETSC_COMM_WORLD, &u);
   VecSetSizes(u, d_dim, PETSC_DECIDE);
   VecSetType(u, VECMPI);
   int offset = d_rank*d_dim;
   int* vec_locs = new int [d_dim];
   for (int i = 0; i < d_dim; ++i) {
      vec_locs[i] = offset + i;
   }
   VecSetValues(u, d_dim, vec_locs, u_in, INSERT_VALUES);
   VecAssemblyBegin(u);
   VecAssemblyEnd(u);
   delete [] vec_locs;

   // If this is the first SVD then build it.  Otherwise add this increment to
   // the system.
   if (d_num_increments == 0) {
      buildInitialSVD(u);
   }
   else {
      buildIncrementalSVD(u);
   }

#ifdef DEBUG
   MatView(d_S, PETSC_VIEWER_STDOUT_WORLD);
   MatView(d_L, PETSC_VIEWER_STDOUT_WORLD);
   MatView(d_U, PETSC_VIEWER_STDOUT_WORLD);
#endif

   // Clean up.
   VecDestroy(&u);
}

void
incremental_svd::buildInitialSVD(
   Vec u)
{
   double one = 1.0;
   int zero = 0;

   // Build d_S.
   MatCreate(PETSC_COMM_WORLD, &d_S);
   MatSetSizes(d_S, PETSC_DECIDE, PETSC_DECIDE, 1, 1);
   MatSetType(d_S, MATDENSE);
   MatSetUp(d_S);
   PetscScalar norm_u;
   VecDot(u, u, &norm_u);
   norm_u = sqrt(norm_u);
   if (d_rank == 0) {
      MatSetValues(d_S, 1, &zero, 1, &zero, &norm_u, INSERT_VALUES);
   }
   MatAssemblyBegin(d_S, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(d_S, MAT_FINAL_ASSEMBLY);

   // Build d_L.
   MatCreate(PETSC_COMM_WORLD, &d_L);
   MatSetSizes(d_L, PETSC_DECIDE, PETSC_DECIDE, 1, 1);
//   MatSetType(d_L, MATDENSE);
   MatSetUp(d_L);
   if (d_rank == 0) {
      MatSetValues(d_L, 1, &zero, 1, &zero, &one, INSERT_VALUES);
   }
   MatAssemblyBegin(d_L, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(d_L, MAT_FINAL_ASSEMBLY);

   // Build d_U.
   MatCreate(PETSC_COMM_WORLD, &d_U);
   MatSetSizes(d_U, d_dim, PETSC_DECIDE, PETSC_DETERMINE, 1);
//   MatSetType(d_U, MATDENSE);
   MatSetUp(d_U);
   int vec_start, vec_end;
   VecGetOwnershipRange(u, &vec_start, &vec_end);
   for (int row = vec_start; row < vec_end; ++row) {
     double u_val;
     VecGetValues(u, 1, &row, &u_val);
     u_val = u_val/norm_u;
     MatSetValues(d_U, 1, &row, 1, &zero, &u_val, INSERT_VALUES);
   }
   MatAssemblyBegin(d_U, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(d_U, MAT_FINAL_ASSEMBLY);

   // We now have another increment.
   ++d_num_increments;
}

void
incremental_svd::buildIncrementalSVD(
   Vec u)
{
   // Compute j, P and the norm of J.
   Vec j;
   Mat P;
   double norm_j;
   compute_J_P_normJ(u, j, P, norm_j);

   // l = P' * u
   Vec l;
   VecCreate(PETSC_COMM_WORLD, &l);
   VecSetSizes(l, PETSC_DECIDE, d_num_increments);
   VecSetType(l, VECMPI);
   MatMultTranspose(P, u, l);

   // Compute cm = u.u and see if this increment is new.
   double cm;
   VecDot(u, u, &cm);

   double eps_squared = d_epsilon*d_epsilon;
   bool is_new_increment = ((cm > eps_squared) &&
      ((norm_j*norm_j)/(eps_squared+cm)) > eps_squared);

   // If this increment is not new and we are skipping redundant increments
   // then clean up and return.
   if (!is_new_increment) {
      if (d_skip_redundant) {
         VecDestroy(&j);
         VecDestroy(&l);
         return;
      }
      else {
         norm_j = 0.0;
      }
   }

   // On process 0, create Q.
   Mat Q;
   constructQ(Q, l, norm_j);

   // Done with l.
   VecDestroy(&l);

   // Now get the singular value decomposition of Q.
   Mat A, sigma, B;
   svd(Q, A, sigma, B);

   // Done with Q so destroy it on process 0 where is was created.
   if (d_rank == 0) {
      MatDestroy(&Q);
   }

   // At this point we either have a new or a redundant increment and we are
   // not skipping redundant increments.
   if (!is_new_increment && !d_skip_redundant) {
      // This increment is not new and we are not skipping redundant
      // increments.
      addRedundantIncrement(A, sigma);
   }
   else if (is_new_increment) {
      // This increment is new.

      // Get the rows of P owned by this process.
      int p_row_start, p_row_end;
      MatGetOwnershipRange(P, &p_row_start, &p_row_end);

      addNewIncrement(j, A, sigma, p_row_start, p_row_end);
   }

   // Clean up.
   MatDestroy(&P);
   VecDestroy(&j);
   MatDestroy(&A);
   MatDestroy(&sigma);
   MatDestroy(&B);
}

void
incremental_svd::compute_J_P_normJ(
   Vec u,
   Vec& j,
   Mat& P,
   double& norm_j)
{
   // j = u
   VecCreate(PETSC_COMM_WORLD, &j);
   VecSetSizes(j, d_dim, PETSC_DECIDE);
   VecSetType(j, VECMPI);
   VecCopy(u, j);

   // P = d_U * d_L
   MatMatMult(d_U, d_L, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &P);

   // Use modified Gram-Schmidt orthogonalization to modify j.
   orthogonalizeJAndComputeNorm(j, P, norm_j);
}

void
incremental_svd::orthogonalizeJAndComputeNorm(
   Vec j,
   Mat P,
   double& norm_j)
{
   // Get the rows of P owned by this process.
   int p_row_start, p_row_end;
   MatGetOwnershipRange(P, &p_row_start, &p_row_end);

   Vec Pvec;
   VecCreate(PETSC_COMM_WORLD, &Pvec);
   VecSetSizes(Pvec, d_dim, PETSC_DECIDE);
   VecSetType(Pvec, VECMPI);
   for (int col = 0; col < d_num_increments; ++col) {
      double factor = 0.0;
      MatGetColumnVector(P, Pvec, col);
      VecDot(Pvec, j, &factor);
      for (int row = p_row_start; row < p_row_end; ++row) {
         double pval, jval;
         VecGetValues(Pvec, 1, &row, &pval);
         VecGetValues(j, 1, &row, &jval);
         double new_jval = jval - factor * pval;
         VecSetValue(j, row, new_jval, INSERT_VALUES);
      }
   }
   VecDestroy(&Pvec);

   // norm_j = sqrt(j.j)
   VecDot(j, j, &norm_j);
   norm_j = sqrt(norm_j);

   // Divide each value of j by norm_j.
   for (int row = p_row_start; row < p_row_end; ++row) {
      double jval;
      VecGetValues(j, 1, &row, &jval);
      jval /= norm_j;
      VecSetValue(j, row, jval, INSERT_VALUES);
   }
}

void
incremental_svd::constructQ(
   Mat& Q,
   Vec l,
   double norm_j)
{
   // Get the rows of d_S owned by this process.
   int d_S_row_start, d_S_row_end;
   MatGetOwnershipRange(d_S, &d_S_row_start, &d_S_row_end);

   // Create an array holding the global ids of all rows in both d_S and Q.
   int* cols = new int [d_num_increments+1];
   for (int i = 0; i < d_num_increments+1; ++i) {
      cols[i] = i;
   }

   // Process 0 is the only process which actually constructs Q.  Q is a small
   // matrix and we only want to get the singular value decomposition of it.
   // Given its size, this can easily and likely more efficiently done on one
   // processor.
   // However, other processors may hold data needed to construct Q so they all
   // must participate in the constuction sending the portions of Q that they
   // hold to processor 0.
   if (d_rank == 0) {
      // Create Q.
      MatCreate(PETSC_COMM_SELF, &Q);
      MatSetSizes(Q, PETSC_DECIDE, PETSC_DECIDE,
                  d_num_increments+1, d_num_increments+1);
      MatSetType(Q, MATDENSE);
      MatSetUp(Q);

      // Add to Q the information held by procesor 0.
      double* qvals = new double [d_num_increments+1];
      for (int row = d_S_row_start; row < d_S_row_end; ++row) {
         MatGetValues(d_S, 1, &row, d_num_increments, cols, qvals);
         VecGetValues(l, 1, &row, &qvals[d_num_increments]);
         MatSetValues(Q, 1, &row, d_num_increments+1, cols,
                      qvals, INSERT_VALUES);
      }
      for (int i = 0; i < d_num_increments; ++i) {
         qvals[i] = 0.0;
      }
      qvals[d_num_increments] = norm_j;
      MatSetValues(Q, 1, &d_num_increments, d_num_increments+1, cols,
                   qvals, INSERT_VALUES);
      delete [] qvals;

      // Get the parts of Q owned by other processors and insert them in Q.
      int size;
      MPI_Comm_size(PETSC_COMM_WORLD, &size);
      const int* ranges;
      MatGetOwnershipRanges(d_S, &ranges);
      MPI_Status status;
      for (int i = 1; i < size; ++i) {
         int num_rows_owned = ranges[i+1] - ranges[i];
         if (num_rows_owned > 0) {
            int num_qvals = (d_num_increments+1) * num_rows_owned;
            qvals = new double [num_qvals];
            int* rows_owned = new int [num_rows_owned];
            for (int row = 0; row < num_rows_owned; ++row) {
               rows_owned[row] = ranges[i] + row;
            }
            MPI_Recv(qvals, num_qvals, MPI_DOUBLE, i, 0,
                     PETSC_COMM_WORLD, &status);
            MatSetValues(Q, num_rows_owned, rows_owned,
                         d_num_increments+1, cols, qvals, INSERT_VALUES);
            delete [] qvals;
            delete [] rows_owned;
         }
      }
   }
   else {
      // Get the parts of Q owned by processors other than 0 and send them to
      // processor 0 for insertion into Q.
      int num_rows_owned = d_S_row_end - d_S_row_start;
      if (num_rows_owned > 0) {
         int num_qvals = (d_num_increments+1) * num_rows_owned;
         double *qvals = new double [num_qvals];
         int offset = 0;
         for (int row = d_S_row_start; row < d_S_row_end; ++row) {
            MatGetValues(d_S, 1, &row, d_num_increments, cols,
                         &qvals[offset*(d_num_increments+1)]);
            VecGetValues(l, 1, &row,
                         &qvals[((offset+1)*d_num_increments)+offset]);
            ++offset;
         }
         MPI_Send(qvals, num_qvals, MPI_DOUBLE, 0, 0, PETSC_COMM_WORLD);
         delete [] qvals;
      }
   }
   delete [] cols;

   // Process 0 assembles Q.
   if (d_rank == 0) {
      MatAssemblyBegin(Q, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(Q, MAT_FINAL_ASSEMBLY);
   }
}

void
incremental_svd::svd(
   Mat A,
   Mat& U,
   Mat& S,
   Mat& V)
{
   // Construct U.
   MatCreate(PETSC_COMM_WORLD, &U);
   MatSetSizes(U, PETSC_DECIDE, PETSC_DECIDE,
               d_num_increments+1, d_num_increments+1);
//   MatSetType(U, MATDENSE);
   MatSetUp(U);

   // Construct S.
   MatCreate(PETSC_COMM_WORLD, &S);
   MatSetSizes(S, PETSC_DECIDE, PETSC_DECIDE,
               d_num_increments+1, d_num_increments+1);
   MatSetType(S, MATDENSE);
   MatSetUp(S);

   // Construct V.
   MatCreate(PETSC_COMM_WORLD, &V);
   MatSetSizes(V, PETSC_DECIDE, PETSC_DECIDE,
               d_num_increments+1, d_num_increments+1);
   MatSetType(V, MATDENSE);
   MatSetUp(V);

   // Q is small and owned only by processor 0 so it is the only process that
   // needs to perform the singular value decomposition.  It will insert the
   // necessary values into U, S, and V.
   if (d_rank == 0) {
      // Solve for the singular value decomposition.
      SVD svd;
      SVDCreate(PETSC_COMM_SELF, &svd);
      SVDSetOperator(svd, A);
      SVDSetType(svd, SVDLAPACK);
      SVDSolve(svd);
      int num_sing_vals;
      SVDGetConverged(svd, &num_sing_vals);
      assert(num_sing_vals == d_num_increments+1);

      // The SLEPc SVD solver does not explicitly export U, S, and V so we
      // need to query the solver for their values and construct them ourself.
      // Note that the solver holds the transpose of V so we need be careful
      // when constructing it.
      Vec u, v;
      MatGetVecs(A, &v, &u);

      int* rows = new int[d_num_increments+1];
      for (int i = 0; i < d_num_increments+1; ++i) {
         rows[i] = i;
      }
      double* uvec = new double [d_num_increments+1];
      double* vvec = new double [d_num_increments+1];
      double* sigma = new double[d_num_increments+1];
      for (int i = 0; i < d_num_increments+1; ++i) {
         sigma[i] = 0.0;
      }
      for (int col = 0; col < d_num_increments+1; ++col) {;
         SVDGetSingularTriplet(svd, col, &sigma[col], u, v);
         VecGetValues(u, d_num_increments+1, rows, uvec);
         VecGetValues(v, d_num_increments+1, rows, vvec);
         MatSetValues(U, d_num_increments+1, rows, 1, &col,
                      uvec, INSERT_VALUES);
         MatSetValues(S, d_num_increments+1, rows, 1, &col,
                      sigma, INSERT_VALUES);
         MatSetValues(V, d_num_increments+1, rows, 1, &col,
                      vvec, INSERT_VALUES);
         sigma[col] = 0.0;
      }
      delete [] rows;
      delete [] uvec;
      delete [] vvec;
      delete [] sigma;
      SVDDestroy(&svd);
      VecDestroy(&u);
      VecDestroy(&v);
   }

   // Assemble U, S, and V.
   MatAssemblyBegin(U, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(U, MAT_FINAL_ASSEMBLY);
   MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);
   MatAssemblyBegin(V, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(V, MAT_FINAL_ASSEMBLY);

   // We want to return the transpose of V, not V.
   MatTranspose(V, MAT_REUSE_MATRIX, &V);
}

void
incremental_svd::addRedundantIncrement(
   Mat A,
   Mat sigma)
{
   // Chop a row and a column off of A to form Amod.  Also form d_S by chopping
   // a row and a column off of sigma.
   int A_row_start, A_row_end;
   MatGetOwnershipRange(A, &A_row_start, &A_row_end);
   Mat Amod;
   MatCreate(PETSC_COMM_WORLD, &Amod);
   MatSetSizes(Amod, PETSC_DECIDE, PETSC_DECIDE,
               d_num_increments, d_num_increments);
//   MatSetType(Amod, MATDENSE);
   MatSetUp(Amod);
   int* cols = new int [d_num_increments];
   for (int i = 0; i < d_num_increments; ++i) {
      cols[i] = i;
   }
   double* vals = new double [d_num_increments];
   for (int row = A_row_start; row < A_row_end; ++row) {
      if (row != d_num_increments) {
         MatGetValues(A, 1, &row, d_num_increments, cols, vals);
         MatSetValues(Amod, 1, &row, d_num_increments, cols,
                      vals, INSERT_VALUES);
         MatGetValues(sigma, 1, &row, d_num_increments, cols, vals);
         MatSetValues(d_S, 1, &row, d_num_increments, cols,
                      vals, INSERT_VALUES);
      }
   }
   MatAssemblyBegin(Amod, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(Amod, MAT_FINAL_ASSEMBLY);
   MatAssemblyBegin(d_S, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(d_S, MAT_FINAL_ASSEMBLY);

   // Multiply d_L and Amod and put result into d_L.
   Mat L_times_Amod;
   MatMatMult(d_L, Amod, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &L_times_Amod);
   MatDestroy(&d_L);
   MatConvert(L_times_Amod, MATSAME, MAT_INITIAL_MATRIX, &d_L);

   // Clean up.
   MatDestroy(&Amod);
   MatDestroy(&L_times_Amod);
   delete [] cols;
   delete [] vals;
}

void
incremental_svd::addNewIncrement(
   Vec j,
   Mat A,
   Mat sigma,
   int d_U_row_start,
   int d_U_row_end)
{
   // Add j as a a new column of d_U.
   Mat newU;
   MatCreate(PETSC_COMM_WORLD, &newU);
   MatSetSizes(newU, d_dim, PETSC_DECIDE, PETSC_DETERMINE, d_num_increments+1);
//   MatSetType(newU, MATDENSE);
   MatSetUp(newU);
   int* cols = new int [d_num_increments+1];
   for (int i = 0; i < d_num_increments+1; ++i) {
      cols[i] = i;
   }
   double* vals = new double [d_num_increments+1];
   for (int row = d_U_row_start; row < d_U_row_end; ++row) {
      MatGetValues(d_U, 1, &row, d_num_increments, cols, vals);
      VecGetValues(j, 1, &row, &vals[d_num_increments]);
      MatSetValues(newU, 1, &row, d_num_increments+1, cols,
                   vals, INSERT_VALUES);
   }
   MatAssemblyBegin(newU, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(newU, MAT_FINAL_ASSEMBLY);
   MatDestroy(&d_U);
   MatConvert(newU, MATSAME, MAT_INITIAL_MATRIX, &d_U);
   MatDestroy(&newU);

   // Add another row and column to d_L.  Only the last value in the new
   // row/column is non-zero and it is 1.
   Mat newL;
   MatCreate(PETSC_COMM_WORLD, &newL);
   MatSetSizes(newL, PETSC_DECIDE, PETSC_DECIDE,
               d_num_increments+1, d_num_increments+1);
//   MatSetType(newL, MATDENSE);
   MatSetUp(newL);
   int d_L_row_start, d_L_row_end;
   MatGetOwnershipRange(d_L, &d_L_row_start, &d_L_row_end);
   for (int row = d_L_row_start; row < d_L_row_end; ++row) {
      MatGetValues(d_L, 1, &row, d_num_increments, cols, vals);
      vals[d_num_increments] = 0.0;
      MatSetValues(newL, 1, &row, d_num_increments+1, cols,
                   vals, INSERT_VALUES);
   }
   if (d_rank == 0) {
      for (int i = 0; i < d_num_increments; ++i) {
         vals[i] = 0.0;
      }
      vals[d_num_increments] = 1.0;
      MatSetValues(newL, 1, &d_num_increments, d_num_increments+1, cols,
                   vals, INSERT_VALUES);
   }
   MatAssemblyBegin(newL, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(newL, MAT_FINAL_ASSEMBLY);
   MatDestroy(&d_L);
   MatMatMult(newL, A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &d_L);
   MatDestroy(&newL);

   // d_S = sigma.
   MatDestroy(&d_S);
   MatConvert(sigma, MATSAME, MAT_INITIAL_MATRIX, &d_S);

   // We now have another increment.
   ++d_num_increments;

   // Clean up.
   delete [] cols;
   delete [] vals;
}

double
incremental_svd::computeNormJ(
   Vec u)
{
   Vec j;
   Mat P;
   double norm_j;
   compute_J_P_normJ(u, j, P, norm_j);
   VecDestroy(&j);
   MatDestroy(&P);
   return norm_j;
}

}
