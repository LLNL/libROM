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

// Description: Interface to the DEIM algorithm to determine the rows of the
// rhs to be sampled for the interpolation of the rhs.

#include "Matrix.h"
#include "mpi.h"
#include <cmath>

namespace CAROM {

// Struct to hold the local maximum absolute value of a basis vector, the row
// it is in, and the processor that owns it.  We will reduce this to find the
// global row containing the maximum of a basis vector.
typedef struct
{
   double row_val;
   int row;
   int proc;
} RowInfo;

// The function to use as an MPI_Op in the reduction to determine the row and
// processor owning the row of the absolute maximum of a basis vector.
void
RowInfoMax(RowInfo* a, RowInfo* b, int* len, MPI_Datatype* type)
{
   for (int i = 0; i < *len; ++i) {
      if (a[i].row_val > b[i].row_val) {
         b[i].row_val = a[i].row_val;
         b[i].row = a[i].row;
         b[i].proc = a[i].proc;
      }
   }
}

void
DEIM(const Matrix* f_basis,
     int num_f_basis_vectors_used,
     int* f_sampled_row,
     int* f_sampled_row_owner,
     Matrix& f_basis_sampled,
     int myid)
{
   // This algorithm determines the rows of f that should be sampled, the
   // processor that owns each sampled row, and fills f_basis_sampled with the
   // sampled rows of the basis of the RHS.

   // Create an MPI_Datatype for the RowInfo struct.
   MPI_Datatype MaxRowType, oldtypes[2];
   int blockcounts[2];
   MPI_Aint offsets[2], extent;
   MPI_Status stat;
   offsets[0] = 0;
   oldtypes[0] = MPI_DOUBLE;
   blockcounts[0] = 1;
   MPI_Type_extent(MPI_DOUBLE, &extent);
   offsets[1] = extent;
   oldtypes[1] = MPI_INT;
   blockcounts[1] = 2;
   MPI_Type_struct(2, blockcounts, offsets, oldtypes, &MaxRowType);
   MPI_Type_commit(&MaxRowType);

   // Create an MPI_Op for the RowInfoMax function.
   MPI_Op RowInfoOp;
   MPI_Op_create((MPI_User_function*)RowInfoMax, true, &RowInfoOp);

   // Get the number of basis vectors and the size of each basis vector.
   int num_basis_vectors = std::min(num_f_basis_vectors_used, f_basis->numColumns());
   int basis_size = f_basis->numRows();

   // The small matrix inverted by the algorithm.  We'll allocate the largest
   // matrix we'll need and set its size at each step in the algorithm.
   Matrix M(num_basis_vectors, num_basis_vectors, false);

   // Scratch space used throughout the algorithm.
   double* c = new double [num_basis_vectors];

   // Figure out the 1st sampled row of the RHS.
   RowInfo f_bv_max_local, f_bv_max_global;
   f_bv_max_local.row_val = -1.0;
   f_bv_max_local.proc = myid;
   for (int i = 0; i < basis_size; ++i) {
      double f_bv_val = fabs(f_basis->item(i, 0));
      if (f_bv_val > f_bv_max_local.row_val) {
         f_bv_max_local.row_val = f_bv_val;
         f_bv_max_local.row = i;
      }
   }
   MPI_Allreduce(&f_bv_max_local, &f_bv_max_global, 1,
                 MaxRowType, RowInfoOp, MPI_COMM_WORLD);
   f_sampled_row[0] = f_bv_max_global.row;
   f_sampled_row_owner[0] = f_bv_max_global.proc;

   // Now get the first sampled row of the basis of the RHS.
   if (f_sampled_row_owner[0] == myid) {
      for (int j = 0; j < num_basis_vectors; ++j) {
         c[j] = f_basis->item(f_sampled_row[0], j);
      }
      MPI_Bcast(c, num_basis_vectors, MPI_DOUBLE,
                f_sampled_row_owner[0], MPI_COMM_WORLD);
   }
   else {
      MPI_Bcast(c, num_basis_vectors, MPI_DOUBLE,
                f_sampled_row_owner[0], MPI_COMM_WORLD);
   }
   // Now add the first sampled row of the basis of the RHS to f_basis_sampled.
   for (int j = 0; j < num_basis_vectors; ++j) {
      f_basis_sampled.item(0, j) = c[j];
   }

   // Now repeat the process for the other sampled rows of the basis of the
   // RHS.
   for (int i = 1; i < num_basis_vectors; ++i) {
      // If we currently know about S sampled rows of the basis of the RHS then
      // M contains the first S columns of those S sampled rows.
      M.setSize(i, i);
      for (int row = 0; row < i; ++row) {
         for (int col = 0; col < i; ++col) {
            M.item(row, col) = f_basis_sampled.item(row, col);
         }
      }

      // Invert M.
      M.inverse();

      // Now compute c, the inverse of M times the next column of the sampled
      // rows of the basis of the RHS.
      for (int minv_row = 0; minv_row < i; ++minv_row) {
         double tmp = 0.0;
         for (int minv_col = 0; minv_col < i; ++minv_col) {
            tmp += M.item(minv_row, minv_col)*
                   f_basis_sampled.item(minv_col, i);
         }
         c[minv_row] = tmp;
      }

      // Now figure out the next sampled row of the basis of f.
      // Compute the first S basis vectors of the RHS times c and find the
      // row of this product have the greatest absolute value.  This is the
      // next sampled row of the basis of f.
      f_bv_max_local.row_val = -1.0;
      f_bv_max_local.proc = myid;
      for (int F_row = 0; F_row < basis_size; ++F_row) {
         double tmp = 0.0;
         for (int F_col = 0; F_col < i; ++F_col) {
            tmp += f_basis->item(F_row, F_col)*c[F_col];
         }
         double r_val = fabs(f_basis->item(F_row, i) - tmp);
         if (r_val > f_bv_max_local.row_val) {
            f_bv_max_local.row_val = r_val;
            f_bv_max_local.row = F_row;
         }
      }
      MPI_Allreduce(&f_bv_max_local, &f_bv_max_global, 1,
                    MaxRowType, RowInfoOp, MPI_COMM_WORLD);
      f_sampled_row[i] = f_bv_max_global.row;
      f_sampled_row_owner[i] = f_bv_max_global.proc;

      // Now get the next sampled row of the basis of f.
      if (f_sampled_row_owner[i] == myid) {
         for (int j = 0; j < num_basis_vectors; ++j) {
            c[j] = f_basis->item(f_sampled_row[i], j);
         }
         MPI_Bcast(c, num_basis_vectors, MPI_DOUBLE,
                   f_sampled_row_owner[i], MPI_COMM_WORLD);
      }
      else {
         MPI_Bcast(c, num_basis_vectors, MPI_DOUBLE,
                   f_sampled_row_owner[i], MPI_COMM_WORLD);
      }
      // Now add the ith sampled row of the basis of the RHS to
      // f_basis_sampled.
      for (int j = 0; j < num_basis_vectors; ++j) {
         f_basis_sampled.item(i, j) = c[j];
      }
   }

   // Free the MPI_Datatype and MPI_Op.
   MPI_Type_free(&MaxRowType);
   MPI_Op_free(&RowInfoOp);

   delete [] c;
}

}
