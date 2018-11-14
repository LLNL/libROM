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

// Description: A class implementing interface of SVD for the static SVD
//              algorithm.

#include "StaticSVD.h"

#include "mpi.h"

#include <stdio.h>
#include <string.h>

/* Use Autotools-detected Fortran name-mangling scheme */
#define dgesdd FC_FUNC(dgesdd, DGESDD)

extern "C" {
void dgesdd(char*, int*, int*, double*, int*,
             double*, double*, int*, double*, int*,
             double*, int*, int*, int*);
}

namespace CAROM {

const int StaticSVD::COMMUNICATE_A = 999;

StaticSVD::StaticSVD(
   int dim,
   int samples_per_time_interval,
   bool debug_algorithm) :
   SVD(dim, samples_per_time_interval, debug_algorithm),
   d_samples(0),
   d_V(0),
   d_this_interval_basis_current(false)
{
   CAROM_ASSERT(dim > 0);
   CAROM_ASSERT(samples_per_time_interval > 0);

   // Get the rank of this process, and the number of processors.
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init) {
      MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
      MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
   }
   else {
      d_rank = 0;
      d_num_procs = 1;
   }

}

StaticSVD::~StaticSVD()
{
   // Delete data members.
   if (d_V) {
      delete d_V;
   }
   for (int i = 0; i < static_cast<int>(d_samples.size()); ++i) {
      if (d_samples[i]) {
         delete [] d_samples[i];
      }
   }
}

bool
StaticSVD::takeSample(
   double* u_in,
   double time,
   bool add_without_increase)
{
   CAROM_ASSERT(u_in != 0);
   CAROM_ASSERT(time >= 0.0);

   // Check the u_in is not non-zero.
   Vector u_vec(u_in, d_dim, true);
   if (u_vec.norm() == 0.0) {
      return false;
   }

   if (isNewTimeInterval()) {

      // We have a new time interval.
      int num_time_intervals =
         static_cast<int>(d_time_interval_start_times.size());
      if (num_time_intervals > 0) {
         if (d_basis) {
            delete d_basis;
            d_basis = 0;
         }
         if (d_basis_right) {
            delete d_basis_right;
            d_basis_right = 0;
         }
         for (int i = 0; i < static_cast<int>(d_samples.size()); ++i) {
            if (d_samples[i]) {
               delete [] d_samples[i];
            }
         }
         d_samples.resize(0);
         delete d_U;
         d_U = 0;
         delete d_S;
         d_S = 0;
         delete d_V;
         d_V = 0;
      }
      d_num_samples = 0;
      d_time_interval_start_times.resize(num_time_intervals+1);
      d_time_interval_start_times[num_time_intervals] = time;
      d_basis = 0;
      d_basis_right = 0;
   }
   double* sample = new double [d_dim];
   memcpy(sample, u_in, d_dim*sizeof(double));
   d_samples.push_back(sample);
   ++d_num_samples;
   d_this_interval_basis_current = false;
   return true;
}

const Matrix*
StaticSVD::getSpatialBasis()
{
   // If this basis is for the last time interval then it may not be up to date
   // so recompute it.
   if (!thisIntervalBasisCurrent()) {
      if (d_basis != 0) {
         delete d_basis;
      }
      if (d_basis_right != 0) {
         delete d_basis_right;
      }
      if (d_U != 0) {
         delete d_U;
      }
      if (d_S != 0) {
         delete d_S;
      }
      if (d_V != 0) {
         delete d_V;
      }
      computeSVD();
   }
   else {
      CAROM_ASSERT(d_basis != 0);
   }
   CAROM_ASSERT(thisIntervalBasisCurrent());
   return d_basis;
}

const Matrix*
StaticSVD::getTemporalBasis()
{
   // If this basis is for the last time interval then it may not be up to date
   // so recompute it.
   if (!thisIntervalBasisCurrent()) {
      if (d_basis != 0) {
         delete d_basis;
      }
      if (d_basis_right != 0) {
         delete d_basis_right;
      }
      if (d_U != 0) {
         delete d_U;
      }
      if (d_S != 0) {
         delete d_S;
      }
      if (d_V != 0) {
         delete d_V;
      }
      computeSVD();
   }
   else {
      CAROM_ASSERT(d_basis_right != 0);
   }
   CAROM_ASSERT(thisIntervalBasisCurrent());
   return d_basis_right;
}

const Matrix*
StaticSVD::getSingularValues()
{
   // If these singular values are for the last time interval then they may not
   // be up to date so recompute them.
   if (!thisIntervalBasisCurrent()) {
      if (d_basis != 0) {
         delete d_basis;
      }
      if (d_basis_right != 0) {
         delete d_basis_right;
      }
      if (d_U != 0) {
         delete d_U;
      }
      if (d_S != 0) {
         delete d_S;
      }
      if (d_V != 0) {
         delete d_V;
      }
      computeSVD();
   }
   else {
      CAROM_ASSERT(d_S != 0);
   }
   CAROM_ASSERT(thisIntervalBasisCurrent());
   return d_S;
}

void
StaticSVD::computeSVD()
{
   // Get the dimensions from each process and the total dimension.
   int* dims = new int [d_num_procs];
   if (d_num_procs > 1) {
      MPI_Allgather(&d_dim, 1, MPI_INT, dims, 1, MPI_INT, MPI_COMM_WORLD);
   }
   else {
      dims[0] = d_dim;
   }

   int total_dim = 0;
   for (int i = 0; i < d_num_procs; ++i) {
      total_dim += dims[i];
   }

   // Do this computation on process 0 and broadcast results to the other
   // processes.
   int num_cols = static_cast<int>(d_samples.size());
   if (d_rank == 0) {
      // Construct storage for the globalized A.
      double* A = new double [total_dim*num_cols];

      // Put this processor's contribution to A into it in column major order.
      int idx = 0;
      for (int col = 0; col < num_cols; ++col) {
         double* col_vals = d_samples[col];
         for (int row = 0; row < d_dim; ++row) {
            A[idx+row] = col_vals[row];
         }
         idx += total_dim;
      }

      // Get the contributions to the globalized A from the other processes.
      if (d_num_procs > 1) {
         int offset = dims[0];
         for (int proc = 1; proc < d_num_procs; ++proc) {
            int this_proc_dim = dims[proc];
            double* Aproc = new double [this_proc_dim*num_cols];
            MPI_Status status;
            MPI_Recv(Aproc,
               this_proc_dim*num_cols,
               MPI_DOUBLE,
               proc,
               COMMUNICATE_A,
               MPI_COMM_WORLD,
               &status);
            int Aproc_idx = 0;
            int Aidx = offset;
            for (int col = 0; col < num_cols; ++col) {
               for (int row = 0; row < this_proc_dim; ++row) {
                  A[Aidx+row] = Aproc[Aproc_idx++];
               }
               Aidx += total_dim;
            }
            delete [] Aproc;
            offset += this_proc_dim;
         }
      }

      // Perform the svd.
      svd(A, total_dim);

      // Broadcast the results.
      if (d_num_procs > 1) {
         MPI_Bcast(&d_U->item(0, 0),
            total_dim*num_cols,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD);
         MPI_Bcast(&d_S->item(0, 0),
            num_cols*num_cols,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD);
         MPI_Bcast(&d_V->item(0, 0),
            num_cols*num_cols,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD);
      }

      // Clean up.
      delete [] A;
   }
   else {
      // Put this processor's contribution to the globalized of A into A in
      // column major order.
      double* A = new double [d_dim*num_cols];
      int idx = 0;
      for (int col = 0; col < num_cols; ++col) {
         double* col_vals = d_samples[col];
         for (int row = 0; row < d_dim; ++row) {
            A[idx++] = col_vals[row];
         }
      }

      // Send the contribution to the globalized A to process 0.
      MPI_Request request;
      MPI_Isend(A,
         d_dim*num_cols,
         MPI_DOUBLE,
         0,
         COMMUNICATE_A,
         MPI_COMM_WORLD,
         &request);

      // Allocate d_U, d_S, and d_V.
      d_U = new Matrix(total_dim, num_cols, false);
      d_S = new Matrix(num_cols, num_cols, false);
      d_V = new Matrix(num_cols, num_cols, false);

      // Get the results from process 0.
      MPI_Bcast(&d_U->item(0, 0),
         total_dim*num_cols,
         MPI_DOUBLE,
         0,
         MPI_COMM_WORLD);
      MPI_Bcast(&d_S->item(0, 0),
         num_cols*num_cols,
         MPI_DOUBLE,
         0,
         MPI_COMM_WORLD);
      MPI_Bcast(&d_V->item(0, 0),
         num_cols*num_cols,
         MPI_DOUBLE,
         0,
         MPI_COMM_WORLD);

      // Clean up.
      delete [] A;
   }
   d_basis = new Matrix(*d_U);
   d_basis_right = new Matrix(*d_V);
   d_this_interval_basis_current = true;
   if (d_debug_algorithm && d_rank == 0) {
      for (int row = 0; row < num_cols; ++row) {
         for (int col = 0; col < num_cols; ++col) {
            printf("%.16e ", d_S->item(row, col));
         }
         printf("\n");
      }
      printf("\n");
      for (int row = 0; row < num_cols; ++row) {
         for (int col = 0; col < num_cols; ++col) {
            printf("%.16e ", d_V->item(row, col));
         }
         printf("\n");
      }
      printf("\n");
      for (int row = 0; row < total_dim; ++row) {
         for (int col = 0; col < num_cols; ++col) {
            printf("%.16e ", d_U->item(row, col));
         }
         printf("\n");
      }
      printf("============================================================\n");
   }
   delete [] dims;
}

void
StaticSVD::svd(
   double* A,
   int total_dim)
{
   CAROM_ASSERT(A != 0);
   CAROM_ASSERT(total_dim > 0);

   int num_samples = static_cast<int>(d_samples.size());

   // Construct d_U.
   d_U = new Matrix(total_dim, num_samples, false);

   // Construct d_S.
   d_S = new Matrix(num_samples, num_samples, false);
   for (int row = 0; row < num_samples; ++row) {
      for (int col = 0; col < num_samples; ++col) {
         d_S->item(row, col) = 0.0;
      }
   }

   // Construct d_V.
   d_V = new Matrix(num_samples, num_samples, false);

   // Use lapack's dgesdd Fortran function to perform the svd.  As this is
   // Fortran A and all the computed matrices are in column major order.
   char jobz = 'A';
   int m = total_dim;
   int n = num_samples;
   int lda = m;
   double* sigma = new double [num_samples];
   double* U = new double[m*m];
   int ldu = m;
   int ldv = num_samples;
   int lwork = n*(4*n + 6)+m;
   double* work = new double [lwork];
   int* iwork = new int [8*n];
   int info;
   dgesdd(&jobz,
	  &m,
	  &n,
	  A,
	  &lda,
	  sigma,
	  U,
	  &ldu,
	  &d_V->item(0, 0),
	  &ldv,
	  work,
	  &lwork,
	  iwork,
	  &info);
   CAROM_ASSERT(info == 0);
   delete [] work;
   delete [] iwork;

   // Place sigma into d_S.
   for (int i = 0; i < num_samples; ++i) {
      d_S->item(i, i) = sigma[i];
   }
   delete [] sigma;

   // Take the first n rows of U and place into d_U.  U is in column major
   // order so convert is to row major order as this is done.
   int uidx = 0;
   for (int row = 0; row < n; ++row) {
      for (int col = 0; col < m; ++col) {
         d_U->item(col, row) = U[uidx++];
      }
   }
   delete [] U;

   // d_V is in column major order.  Convert it to row major order.
   /*
   for (int row = 0; row < num_samples; ++row) {
      for (int col = row+1; col < num_samples; ++col) {
         double tmp = d_V->item(row, col);
         d_V->item(row, col) = d_V->item(col, row);
         d_V->item(col, row) = tmp;
      }
   }*/
}

}
