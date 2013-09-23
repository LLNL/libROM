#include "incremental_svd.h"
#include "SLEPcManager.h"
#include "slepc.h"
#include <cmath>
#include <assert.h>

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
   d_S(0)
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
   VecNorm(u, NORM_2, &norm_u);
   if (d_rank == 0) {
      MatSetValues(d_S, 1, &zero, 1, &zero, &norm_u, INSERT_VALUES);
   }
   MatAssemblyBegin(d_S, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(d_S, MAT_FINAL_ASSEMBLY);

   // Build d_L.
   MatCreate(PETSC_COMM_WORLD, &d_L);
   MatSetSizes(d_L, PETSC_DECIDE, PETSC_DECIDE, 1, 1);
   MatSetUp(d_L);
   if (d_rank == 0) {
      MatSetValues(d_L, 1, &zero, 1, &zero, &one, INSERT_VALUES);
   }
   MatAssemblyBegin(d_L, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(d_L, MAT_FINAL_ASSEMBLY);

   // Build d_U.
   MatCreate(PETSC_COMM_WORLD, &d_U);
   MatSetSizes(d_U, d_dim, PETSC_DECIDE, PETSC_DETERMINE, 1);
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
   d_U_Times_d_L(P);

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
   int num_rows = p_row_end - p_row_start;
   int* rows = new int [num_rows];
   for (int i = 0; i < num_rows; ++i) {
      rows[i] = i + p_row_start;
   }
   double* pvec_vals = new double [num_rows];

   Vec Pvec;
   VecCreate(PETSC_COMM_WORLD, &Pvec);
   VecSetSizes(Pvec, d_dim, PETSC_DECIDE);
   VecSetType(Pvec, VECMPI);
   for (int col = 0; col < d_num_increments; ++col) {
      double factor;
      MatGetValues(P, num_rows, rows, 1, &col, pvec_vals);
      VecSetValues(Pvec, num_rows, rows, pvec_vals, INSERT_VALUES);
      VecDot(Pvec, j, &factor);
      VecAXPY(j, -factor, Pvec);
   }
   VecAssemblyBegin(j);
   VecAssemblyEnd(j);
   delete [] rows;
   delete [] pvec_vals;

   // Done with Pvec.
   VecDestroy(&Pvec);

   // Normalize j.
   VecNormalize(j, &norm_j);
}

void
incremental_svd::constructQ(
   Mat& Q,
   Vec l,
   double norm_j)
{
   // Get the rows of d_S owned by each process.
   const int* ranges;
   MatGetOwnershipRanges(d_S, &ranges);

   // Create an array holding the global ids of all rows and columns in both
   // d_S and Q.
   int* row_cols = new int [d_num_increments+1];
   for (int i = 0; i < d_num_increments+1; ++i) {
      row_cols[i] = i;
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
   }

   // Processor 0 allocates storage for all the values of Q.
   double* all_qvals = 0;
   if (d_rank == 0) {
      all_qvals = new double [(d_num_increments+1)*(d_num_increments+1)];
   }

   // Each processor allocates and fills an array for it's contribution to Q.
   int my_num_rows_owned = ranges[d_rank+1] - ranges[d_rank];
   int num_my_qvals = (d_num_increments+1)*my_num_rows_owned;
   double* my_qvals = 0;
   if (num_my_qvals > 0) {
      my_qvals = new double [num_my_qvals];
   }
   int offset = 0;
   for (int row = ranges[d_rank]; row < ranges[d_rank+1]; ++row) {
      MatGetValues(d_S, 1, &row, d_num_increments, row_cols,
                   &my_qvals[offset]);
      VecGetValues(l, 1, &row, &my_qvals[offset+d_num_increments]);
      offset += d_num_increments + 1;
   }

   // Process 0 need to know how much data it will receive from each process
   // and the displacement in all_qvals for the data from each process.
   int* recv_cts = 0;
   int* displ = 0;
   if (d_rank == 0) {
      recv_cts = new int[d_size];
      displ = new int[d_size];
      offset = 0;
      for (int i = 0; i < d_size; ++i) {
         recv_cts[i] = (ranges[i+1] - ranges[i]) * (d_num_increments+1);
         displ[i] = offset;
         offset += recv_cts[i];
      }
   }

   // Each process sends to process 0 its contribution to Q.
   MPI_Gatherv(my_qvals, num_my_qvals, MPI_DOUBLE, all_qvals,
               recv_cts, displ, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

   // Process 0 sets the values of Q and assembles the matrix.
   if (d_rank == 0) {
      offset = d_num_increments*(d_num_increments+1);
      for (int i = 0; i < d_num_increments; ++i) {
         all_qvals[offset+i] = 0.0;
      }
      all_qvals[offset+d_num_increments] = norm_j;
      MatSetValues(Q, d_num_increments+1, row_cols, d_num_increments+1,
                   row_cols, all_qvals, INSERT_VALUES);
      MatAssemblyBegin(Q, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(Q, MAT_FINAL_ASSEMBLY);
   }

   // Clean up.
   if (d_rank == 0) {
      delete [] all_qvals;
      delete [] recv_cts;
      delete [] displ;
   }
   if (my_qvals) {
      delete [] my_qvals;
   }
   delete [] row_cols;
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

      // Clean up.
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
   MatSetType(newU, MATDENSE);
   MatSetUp(newU);
   int* cols = new int [d_num_increments+1];
   for (int i = 0; i < d_num_increments+1; ++i) {
      cols[i] = i;
   }
   int num_d_U_rows = d_U_row_end - d_U_row_start;
   int* rows = new int [num_d_U_rows];
   for (int i = 0; i < num_d_U_rows; ++i) {
      rows[i] = i+d_U_row_start;
   }
   double* vals = new double [num_d_U_rows*d_num_increments];
   MatGetValues(d_U, num_d_U_rows, rows, d_num_increments, cols, vals);
   MatSetValues(newU, num_d_U_rows, rows, d_num_increments, cols,
                vals, INSERT_VALUES);
   VecGetValues(j, num_d_U_rows, rows, vals);
   MatSetValues(newU, num_d_U_rows, rows, 1, &cols[d_num_increments],
                vals, INSERT_VALUES);
   delete [] rows;
   delete [] vals;
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
   MatSetUp(newL);
   int d_L_row_start, d_L_row_end;
   vals = new double [d_num_increments+1];
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
   delete [] vals;
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
}

double
incremental_svd::computeNormJ(
   Vec u)
{
   Vec j;
   Mat P;
   double norm_j;
   compute_J_P_normJ(u, j, P, norm_j);

   // Clean up and return the norm of j.
   VecDestroy(&j);
   MatDestroy(&P);
   return norm_j;
}

void
incremental_svd::d_U_Times_d_L(
   Mat& P)
{
   // Get the part of d_L that this processor owns.
   const int* ranges;
   MatGetOwnershipRanges(d_L, &ranges);
   int my_num_L_rows = ranges[d_rank+1] - ranges[d_rank];
   int my_num_vals = my_num_L_rows*d_num_increments;
   double* my_vals = 0;
   int* cols = new int [d_num_increments];
   for (int i = 0; i < d_num_increments; ++i) {
      cols[i] = i;
   }
   if (my_num_L_rows) {
      int* my_L_rows = new int [my_num_L_rows];
      for (int i = 0; i < my_num_L_rows; ++i) {
         my_L_rows[i] = i + ranges[d_rank];
      }
      my_vals = new double [my_num_vals];
      MatGetValues(d_L, my_num_L_rows, my_L_rows,
                   d_num_increments, cols, my_vals);
      delete [] my_L_rows;
   }

   // Send this processors part of d_L to all other processors so that all
   // processors own a complete version of d_L.
   int* recvcts = new int [d_size];
   int* offsets = new int [d_size];
   offsets[0] = 0;
   for (int i = 0; i < d_size; ++i) {
      recvcts[i] = d_num_increments*(ranges[i+1]-ranges[i]);
      if (i > 0) {
         offsets[i] = offsets[i-1] + recvcts[i-1];
      }
   }
   double* L = new double [d_num_increments*d_num_increments];
   MPI_Allgatherv(my_vals, my_num_vals, MPI_DOUBLE,
                  L, recvcts, offsets, MPI_DOUBLE, PETSC_COMM_WORLD);
   delete [] recvcts;
   delete [] offsets;
   if (my_vals) {
      delete [] my_vals;
   }

   // Get the part of d_U that this processor owns.
   int d_U_row_start, d_U_row_end;
   MatGetOwnershipRange(d_U, &d_U_row_start, &d_U_row_end);
   int num_d_U_rows = d_U_row_end - d_U_row_start;
   int* rows = new int [num_d_U_rows];
   for (int i = 0; i < num_d_U_rows; ++i) {
      rows[i] = i + d_U_row_start;
   }
   double* U = new double [d_dim*d_num_increments];
   MatGetValues(d_U, d_dim, rows, d_num_increments, cols, U);

   // Create the result.
   MatCreate(PETSC_COMM_WORLD, &P);
   MatSetSizes(P, d_dim, PETSC_DECIDE, PETSC_DETERMINE, d_num_increments);
   MatSetType(P, MATDENSE);
   MatSetUp(P);

   // Construct the product.
   double* P_vals = new double [d_dim * d_num_increments];
   for (int L_col = 0; L_col < d_num_increments; ++L_col) {
      int U_idx = 0;
      int P_idx = L_col;
      for (int U_row = 0; U_row < d_dim; ++U_row) {
         int L_idx = L_col;
         P_vals[P_idx] = 0.0;
         for (int entry = 0; entry < d_num_increments; ++entry) {
            P_vals[P_idx] += U[U_idx]*L[L_idx];
            ++U_idx;
            L_idx += d_num_increments;
         }
         P_idx += d_num_increments;
      }
   }

   // Put the product into the result.
   MatSetValues(P, d_dim, rows, d_num_increments, cols, P_vals, INSERT_VALUES);
   MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY);

   // Clean up.
   delete [] rows;
   delete [] cols;
   delete [] U;
   delete [] L;
   delete [] P_vals;
}

}
