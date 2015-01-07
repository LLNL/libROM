/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The abstract incremental SVD algorithm defines algorithm
 *              interface.
 *
 *****************************************************************************/

#include "IncrementalSVD.h"

#include "mpi.h"

#include <cmath>
#include <limits>
#include <stdio.h>

extern "C" {
void dgesdd_(char*, int*, int*, double*, int*,
             double*, double*, int*, double*, int*,
             double*, int*, int*, int*);
}

namespace CAROM {

const int IncrementalSVD::COMMUNICATE_U = 666;

IncrementalSVD::IncrementalSVD(
   int dim,
   double redundancy_tol,
   bool skip_redundant,
   int samples_per_time_interval,
   bool debug_rom) :
   SVD(dim, samples_per_time_interval, debug_rom),
   d_redundancy_tol(redundancy_tol),
   d_skip_redundant(skip_redundant),
   d_S(0),
   d_total_dim(0)
{
   CAROM_ASSERT(redundancy_tol > 0.0);

   // Get the number of processors, the dimensions for each process, and the
   // total dimension.
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init) {
      MPI_Comm_size(MPI_COMM_WORLD, &d_size);
   }
   else {
      d_size = 1;
   }
   d_proc_dims.reserve(d_size);
   if (mpi_init) {
      MPI_Allgather(&d_dim,
         1,
         MPI_INT,
         &d_proc_dims[0],
         1,
         MPI_INT,
         MPI_COMM_WORLD);
   }
   else {
      d_proc_dims[0] = d_dim;
   }
   for (int i = 0; i < d_size; ++i) {
      d_total_dim += d_proc_dims[i];
   }
}

IncrementalSVD::~IncrementalSVD()
{
   // Delete data members.
   if (d_S) {
      delete d_S;
   }
}

void
IncrementalSVD::takeSample(
   const double* u_in,
   double time)
{
   CAROM_ASSERT(u_in != 0);
   CAROM_ASSERT(time >= 0.0);

   // If this is the first SVD then build it.  Otherwise add this sample to the
   // system.
   if (isNewTimeInterval()) {
      buildInitialSVD(u_in, time);
   }
   else {
      buildIncrementalSVD(u_in);
   }

   if (d_debug_rom) {
      int mpi_init;
      MPI_Initialized(&mpi_init);
      int rank;
      if (mpi_init) {
         MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      }
      else {
         rank = 0;
      }
      const Matrix* basis = getBasis();
      if (rank == 0) {
         // Print d_S.
         for (int row = 0; row < d_num_samples; ++row) {
            for (int col = 0; col < d_num_samples; ++col) {
               printf("%.16e  ", d_S->item(row, col));
            }
            printf("\n");
         }
         printf("\n");

         // Print process 0's part of the basis.
         for (int row = 0; row < d_dim; ++row) {
            for (int col = 0; col < d_num_samples; ++col) {
               printf("%.16e ", basis->item(row, col));
            }
            printf("\n");
         }

         // Gather other processor's parts of the basis and print them.
         for (int proc = 1; proc < d_size; ++proc) {
            double* m = new double[d_proc_dims[proc]*d_num_samples];
            MPI_Status status;
            MPI_Recv(m,
               d_proc_dims[proc]*d_num_samples,
               MPI_DOUBLE,
               proc,
               COMMUNICATE_U,
               MPI_COMM_WORLD,
               &status);
            int idx = 0;
            for (int row = 0; row < d_proc_dims[proc]; ++row) {
               for (int col = 0; col < d_num_samples; ++col) {
                  printf("%.16e ", m[idx++]);
               }
               printf("\n");
            }
            delete [] m;
         }
         printf("============================================================\n");
      }
      else {
         // Send this processor's part of the basis to process 0.
         MPI_Request request;
         MPI_Isend(const_cast<double*>(&basis->item(0, 0)),
            d_dim*d_num_samples,
            MPI_DOUBLE,
            0,
            COMMUNICATE_U,
            MPI_COMM_WORLD,
            &request);
      }
   }
}

const Matrix*
IncrementalSVD::getBasis()
{
   CAROM_ASSERT(d_basis != 0);
   return d_basis;
}

void
IncrementalSVD::buildIncrementalSVD(
   const double* u)
{
   CAROM_ASSERT(u != 0);

   // l = basis' * u
   Vector u_vec(u, d_dim, true);
   Vector* l = d_basis->transposeMult(u_vec);

   // basisl = basis * l
   Vector* basisl = d_basis->mult(l);

   // Compute k = sqrt(u.u - 2.0*l.l + basisl.basisl) which is ||u - basisl||.
   double k = u_vec.inner_product(u_vec) - 2.0*l->inner_product(l) +
      basisl->inner_product(basisl);
   if (k <= 0) {
      k = 0;
   }
   else {
      k = sqrt(k);
   }

   // Use k to see if this sample is new.
   double default_tol;
   if (d_num_samples == 1) {
      default_tol = std::numeric_limits<double>::epsilon();
   }
   else {
      default_tol = d_num_samples*std::numeric_limits<double>::epsilon()*d_S->item(0, 0);
   }
   if (d_redundancy_tol < default_tol) {
      d_redundancy_tol = default_tol;
   }
   bool is_new_sample;
   if (k < d_redundancy_tol) {
      k = 0;
      is_new_sample = false;
   }
   else {
      is_new_sample = true;
   }

   // Create Q.
   double* Q;
   constructQ(Q, l, k);

   // Now get the singular value decomposition of Q.
   Matrix* A;
   Matrix* sigma;
   svd(Q, A, sigma);

   // Done with Q.
   delete [] Q;

   // We need to add the sample if it is new or if it is redundant and we are
   // not skipping redundant samples.
   if (!is_new_sample && !d_skip_redundant) {
      // This sample is redundant and we are not skipping redundant samples.
      addRedundantSample(A, sigma);
      delete sigma;
   }
   else if (is_new_sample) {
      // This sample is new.

      // Compute j
      Vector* j = u_vec.minus(basisl);
      for (int i = 0; i < d_dim; ++i) {
         j->item(i) /= k;
      }

      // addNewSample will assign sigma to d_S hence it should not be deleted
      // upon return.
      addNewSample(j, A, sigma);
      delete j;
   }
   delete basisl;

   // Compute the basis vectors.
   delete d_basis;
   computeBasis();

   // Clean up.
   delete l;
   delete A;
}

void
IncrementalSVD::constructQ(
   double*& Q,
   const Vector* l,
   double k)
{
   CAROM_ASSERT(l != 0);
   CAROM_ASSERT(l->dim() == numSamples());

   // Create Q.
   Q = new double [(d_num_samples+1)*(d_num_samples+1)];

   // Fill Q in column major order.
   int q_idx = 0;
   for (int row = 0; row < d_num_samples; ++row) {
      q_idx = row;
      for (int col = 0; col < d_num_samples; ++col) {
         Q[q_idx] = d_S->item(row, col);
         q_idx += d_num_samples+1;
      }
      Q[q_idx] = l->item(row);
   }
   q_idx = d_num_samples;
   for (int col = 0; col < d_num_samples; ++col) {
      Q[q_idx] = 0.0;
      q_idx += d_num_samples+1;
   }
   Q[q_idx] = k;
}

void
IncrementalSVD::svd(
   double* A,
   Matrix*& U,
   Matrix*& S)
{
   CAROM_ASSERT(A != 0);

   // Construct U, S, and V.
   U = new Matrix(d_num_samples+1, d_num_samples+1, false);
   S = new Matrix(d_num_samples+1, d_num_samples+1, false);
   Matrix* V = new Matrix(d_num_samples+1, d_num_samples+1, false);
   for (int row = 0; row < d_num_samples+1; ++row) {
      for (int col = 0; col < d_num_samples+1; ++col) {
         S->item(row, col) = 0.0;
      }
   }

   // Use lapack's dgesdd_ Fortran function to perform the svd.  As this is
   // Fortran A and all the computed matrices are in column major order.
   double* sigma = new double [d_num_samples+1];
   char jobz = 'A';
   int m = d_num_samples+1;
   int n = d_num_samples+1;
   int lda = d_num_samples+1;
   int ldu = d_num_samples+1;
   int ldv = d_num_samples+1;
   int lwork = m*(4*m + 7);
   double* work = new double [lwork];
   int iwork[8*m];
   int info;
   dgesdd_(&jobz,
      &m,
      &n,
      A,
      &lda,
      sigma,
      &U->item(0, 0),
      &ldu,
      &V->item(0, 0),
      &ldv,
      work,
      &lwork,
      iwork,
      &info);
   CAROM_ASSERT(info == 0);
   delete [] work;

   // Place sigma into S.
   for (int i = 0; i < d_num_samples+1; ++i) {
      S->item(i, i) = sigma[i];
   }
   delete [] sigma;

   // U is column major order so convert it to row major order.
   for (int row = 0; row < d_num_samples+1; ++row) {
      for (int col = row+1; col < d_num_samples+1; ++col) {
         double tmp = U->item(row, col);
         U->item(row, col) = U->item(col, row);
         U->item(col, row) = tmp;
      }
   }
   delete V;
}

double
IncrementalSVD::checkOrthogonality(
   const Matrix* m)
{
   CAROM_ASSERT(m != 0);

   double result = 0.0;
   if (d_num_samples > 1) {
      int last_col = d_num_samples-1;
      double tmp = 0.0;
      int num_rows = m->numRows();
      for (int i = 0; i < num_rows; ++i) {
         tmp += m->item(i, 0) * m->item(i, last_col);
      }
      if (m->distributed() && d_size > 1) {
         MPI_Allreduce(&tmp, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
      else {
         result = tmp;
      }
   }
   return result;
}

void
IncrementalSVD::reOrthogonalize(
   Matrix* m)
{
   CAROM_ASSERT(m != 0);

   int num_rows = m->numRows();
   int num_cols = m->numColumns();
   for (int work = 1; work < num_cols; ++work) {
      double tmp;
      for (int col = 0; col < work; ++col) {
         double factor = 0.0;
         tmp = 0.0;
         for (int i = 0; i < num_rows; ++i) {
            tmp += m->item(i, col)*m->item(i, work);
         }
         if (d_size > 1) {
            MPI_Allreduce(&tmp,
               &factor,
               1,
               MPI_DOUBLE,
               MPI_SUM,
               MPI_COMM_WORLD);
         }
         else {
            factor = tmp;
         }

         for (int i = 0; i < num_rows; ++i) {
            m->item(i, work) -= factor*m->item(i, col);
         }
      }
      double norm = 0.0;
      tmp = 0.0;
      for (int i = 0; i < num_rows; ++i) {
         tmp += m->item(i, work)*m->item(i, work);
      }
      if (d_size > 1) {
         MPI_Allreduce(&tmp, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
      else {
         norm = tmp;
      }
      norm = sqrt(norm);
      for (int i = 0; i < num_rows; ++i) {
         m->item(i, work) /= norm;
      }
   }
}

}
