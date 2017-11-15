/******************************************************************************
 *
 * Copyright (c) 2013-2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 * Written by William Arrighi wjarrighi@llnl.gov
 * CODE-686965
 * All rights reserved.
 *
 * This file is part of libROM.
 * For details, see https://computation.llnl.gov/librom
 * Please also read README_BSD_NOTICE.
 *
 * Redistribution and use in source and binary forms, with or without
 * modifications, are permitted provided that the following conditions are met:
 *
 *    o Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the disclaimer below.
 *    o Redistribution in binary form must reproduce the above copyright
 *      notice, this list of conditions and the disclaimer (as noted below) in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *    o Neither the name of the LLNS/LLNL nor the names of its contributors may
 *      be used to endorse or promote products derived from this software
 *      without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
 * LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OR SUCH DAMAGE.
 *
 *****************************************************************************/

// Description: The abstract incremental SVD algorithm defines algorithm
//              interface.

#include "IncrementalSVD.h"
#include "HDFDatabase.h"

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
   double linearity_tol,
   bool skip_linearly_dependent,
   int samples_per_time_interval,
   bool save_state,
   bool restore_state,
   bool debug_algorithm) :
   SVD(dim, samples_per_time_interval, debug_algorithm),
   d_linearity_tol(linearity_tol),
   d_skip_linearly_dependent(skip_linearly_dependent),
   d_total_dim(0),
   d_save_state(save_state),
   d_state_database(0)
{
   CAROM_ASSERT(linearity_tol > 0.0);

   // Get the number of processors, the dimensions for each process, and the
   // total dimension.
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init) {
      MPI_Comm_size(MPI_COMM_WORLD, &d_size);
      MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
   }
   else {
      d_size = 1;
      d_rank = 0;
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

   // If the state of the SVD is to be restored then open the database and
   // restore the necessary data from the database now.
   if (restore_state) {
      // Open state database file.
      char file_name[100];
      sprintf(file_name, "state.%06d", d_rank);
      d_state_database = new HDFDatabase();
      bool is_good = d_state_database->open(file_name);
      if (is_good) {
         // Read time interval start time.
         double time;
         d_state_database->getDouble("time", time);
         d_time_interval_start_times.resize(1);
         d_time_interval_start_times[0] = time;

         // Read d_U.
         int num_rows;
         d_state_database->getInteger("U_num_rows", num_rows);
         int num_cols;
         d_state_database->getInteger("U_num_cols", num_cols);
         d_U = new Matrix(num_rows, num_cols, true);
         d_state_database->getDoubleArray("U",
                                          &d_U->item(0, 0),
                                          num_rows*num_cols);

         // Read d_S.
         d_state_database->getInteger("S_num_rows", num_rows);
         d_state_database->getInteger("S_num_cols", num_cols);
         d_S = new Matrix(num_rows, num_cols, false);
         d_state_database->getDoubleArray("S",
                                          &d_S->item(0, 0),
                                          num_rows*num_cols);

         // Set d_num_samples.
         d_num_samples = num_cols;
      }
      else {
         delete d_state_database;
         d_state_database = 0;
      }
   }
}

IncrementalSVD::~IncrementalSVD()
{
   // If the state of the SVD is to be saved, then save d_S and d_U now.  The
   // derived class has already created the database.
   //
   // If there are multiple time intervals then saving and restoring the state
   // does not make sense as there is not one, all encompassing, basis.
   if (d_save_state && d_time_interval_start_times.size() == 1) {
      // Save the time interval start time.
      d_state_database->putDouble("time", d_time_interval_start_times[0]);

      // Save d_U.
      int num_rows = d_U->numRows();
      d_state_database->putInteger("U_num_rows", num_rows);
      int num_cols = d_U->numColumns();
      d_state_database->putInteger("U_num_cols", num_cols);
      d_state_database->putDoubleArray("U", &d_U->item(0, 0), num_rows*num_cols);

      // Save d_S.
      num_rows = d_S->numRows();
      d_state_database->putInteger("S_num_rows", num_rows);
      num_cols = d_S->numColumns();
      d_state_database->putInteger("S_num_cols", num_cols);
      d_state_database->putDoubleArray("S",
                                       &d_S->item(0, 0),
                                       num_rows*num_cols);

      // Close state database file and delete database object.
      d_state_database->close();
      delete d_state_database;
   }
}

bool
IncrementalSVD::takeSample(
   const double* u_in,
   double time)
{
   CAROM_ASSERT(u_in != 0);
   CAROM_ASSERT(time >= 0.0);

   // Check that u_in is not non-zero.
   Vector u_vec(u_in, d_dim, true);
   if (u_vec.norm() == 0.0) {
      return false;
   }

   // If this is the first SVD then build it.  Otherwise add this sample to the
   // system.
   bool result = true;
   if (isNewTimeInterval()) {
      buildInitialSVD(u_in, time);
   }
   else {
      result = buildIncrementalSVD(u_in);
   }

   if (d_debug_algorithm) {
      const Matrix* basis = getBasis();
      if (d_rank == 0) {
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
   return result;
}

const Matrix*
IncrementalSVD::getBasis()
{
   CAROM_ASSERT(d_basis != 0);
   return d_basis;
}

const Matrix*
IncrementalSVD::getSingularValues()
{
   CAROM_ASSERT(d_S != 0);
   return d_S;
}

bool
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
   // This is the error in the projection of u into the reduced order space
   // and subsequent lifting back to the full order space.
   double k = u_vec.inner_product(u_vec) - 2.0*l->inner_product(l) +
      basisl->inner_product(basisl);
   if (k <= 0) {
      k = 0;
   }
   else {
      k = sqrt(k);
   }

   // Use k to see if this sample is new.
   bool linearly_dependent_sample;
   if (k < d_linearity_tol) {
      k = 0;
      linearly_dependent_sample = true;
   }
   else {
      linearly_dependent_sample = false;
   }

   // Create Q.
   double* Q;
   constructQ(Q, l, k);
   delete l;

   // Now get the singular value decomposition of Q.
   Matrix* A;
   Matrix* sigma;
   bool result = svd(Q, A, sigma);

   // Done with Q.
   delete [] Q;

   // If the svd was successful then add the sample.  Otherwise clean up and
   // return.
   if (result) {

      // We need to add the sample if it is not linearly dependent or if it is
      // linearly dependent and we are not skipping linearly dependent samples.
      if (linearly_dependent_sample && !d_skip_linearly_dependent) {
         // This sample is linearly dependent and we are not skipping linearly
         // dependent samples.
         addLinearlyDependentSample(A, sigma);
         delete sigma;
      }
      else if (!linearly_dependent_sample) {
         // This sample is not linearly dependent.

         // Compute j
         Vector* j = u_vec.minus(basisl);
         for (int i = 0; i < d_dim; ++i) {
            j->item(i) /= k;
         }

         // addNewSample will assign sigma to d_S hence it should not be
         // deleted upon return.
         addNewSample(j, A, sigma);
         delete j;
      }
      delete basisl;
      delete A;

      // Compute the basis vectors.
      delete d_basis;
      computeBasis();
   }
   else {
      delete basisl;
      delete A;
      delete sigma;
   }
   return result;
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

bool
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
   delete [] work;

   // If the svd succeeded, fill U and S.  Otherwise clean up and return.
   if (info == 0) {
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
   else {
      delete V;
      delete [] sigma;
   }
   return info == 0;
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
