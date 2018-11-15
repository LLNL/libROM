/******************************************************************************
 *
 * Copyright (c) 2013-2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 * Written by Dylan Matthew Copeland, copeland11@llnl.gov
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
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include "DEIM.h"

namespace CAROM {

void GNAT(const Matrix* f_basis,
	  const int num_f_basis_vectors_used,
	  int* f_sampled_row,
	  int* f_sampled_rows_per_proc,
	  Matrix& f_basis_sampled_inv,
	  const int myid,
	  const int num_procs,
	  const int num_samples_req)
{
  // This algorithm determines the rows of f that should be sampled, the
  // processor that owns each sampled row, and fills f_basis_sampled_inv with
  // the inverse of the sampled rows of the basis of the RHS.

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
  const int num_basis_vectors =
    std::min(num_f_basis_vectors_used, f_basis->numColumns());
  const int basis_size = f_basis->numRows();

  const int num_samples = num_samples_req > 0 ? num_samples_req : num_basis_vectors;
  
  const int ns_mod_nr = num_samples % num_basis_vectors;
  int ns = 0;
   
  CAROM_ASSERT(num_samples >= num_basis_vectors && num_samples <= basis_size);

  // The small matrix inverted by the algorithm.  We'll allocate the largest
  // matrix we'll need and set its size at each step in the algorithm.
  Matrix M(num_basis_vectors, num_basis_vectors, false);

  // Scratch space used throughout the algorithm.
  double* c = new double [num_basis_vectors];

  std::vector<std::set<int> > proc_sampled_f_row(num_procs);
  std::vector<std::map<int, int> > proc_f_row_to_tmp_fs_row(num_procs);
  int num_f_basis_cols = f_basis_sampled_inv.numColumns();
  Matrix tmp_fs(f_basis_sampled_inv.numRows(),
		num_f_basis_cols,
		f_basis_sampled_inv.distributed());

  // Figure out the 1st sampled rows of the RHS.
  RowInfo f_bv_max_local, f_bv_max_global;

  const int ns0 = 0 < ns_mod_nr ? (num_samples / num_basis_vectors) + 1 : num_samples / num_basis_vectors;

  for (int k=0; k<ns0; ++k)
    {
      f_bv_max_local.row_val = -1.0;
      f_bv_max_local.proc = myid;
      for (int i = 0; i < basis_size; ++i) {
	// Check whether this row has already been sampled.
	std::set<int>::const_iterator found = proc_sampled_f_row[myid].find(i);
	if (found == proc_sampled_f_row[myid].end()) // not found
	  {
	    double f_bv_val = fabs(f_basis->item(i, 0));
	    if (f_bv_val > f_bv_max_local.row_val) {
	      f_bv_max_local.row_val = f_bv_val;
	      f_bv_max_local.row = i;
	    }
	  }
      }
      MPI_Allreduce(&f_bv_max_local, &f_bv_max_global, 1,
		    MaxRowType, RowInfoOp, MPI_COMM_WORLD);

      // Now get the first sampled row of the basis of the RHS.
      if (f_bv_max_global.proc == myid) {
	for (int j = 0; j < num_basis_vectors; ++j) {
	  c[j] = f_basis->item(f_bv_max_global.row, j);
	}
      }
      MPI_Bcast(c, num_basis_vectors, MPI_DOUBLE,
		f_bv_max_global.proc, MPI_COMM_WORLD);
      // Now add the first sampled row of the basis of the RHS to tmp_fs.
      for (int j = 0; j < num_basis_vectors; ++j) {
	tmp_fs.item(k, j) = c[j];
      }
      proc_sampled_f_row[f_bv_max_global.proc].insert(f_bv_max_global.row);
      proc_f_row_to_tmp_fs_row[f_bv_max_global.proc][f_bv_max_global.row] = k;
    }
  
  ns += ns0;
  
  // Now repeat the process for the other sampled rows of the basis of the RHS.
  for (int i = 1; i < num_basis_vectors; ++i) {
    const int nsi = i < ns_mod_nr ? (num_samples / num_basis_vectors) + 1 : num_samples / num_basis_vectors;

    // If we currently know about S sampled rows of the basis of the RHS then
    // M contains the first i columns of those S sampled rows.
    M.setSize(ns, i);
    for (int row = 0; row < ns; ++row) {
      for (int col = 0; col < i; ++col) {
	M.item(row, col) = tmp_fs.item(row, col);
      }
    }

    // Compute the pseudo-inverse of M, storing its transpose.
    M.pseudoinverse();

    // Now compute c, the inverse of M times the next column of the sampled
    // rows of the basis of the RHS.
    for (int minv_row = 0; minv_row < i; ++minv_row) {
      double tmp = 0.0;
      for (int minv_col = 0; minv_col < ns; ++minv_col) {
	if (ns == i)
	  tmp += M.item(minv_row, minv_col)*tmp_fs.item(minv_col, i);
	else
	  tmp += M.item(minv_col, minv_row)*tmp_fs.item(minv_col, i);  // Transposing M^+, which is stored as its transpose.
      }
      c[minv_row] = tmp;
    }
      
    for (int k=0; k<nsi; ++k)
      {
	// Now figure out the next sampled row of the basis of f.
	// Compute the first S basis vectors of the RHS times c and find the
	// row of this product have the greatest absolute value.  This is the
	// next sampled row of the basis of f.
	f_bv_max_local.row_val = -1.0;
	f_bv_max_local.proc = myid;
	for (int F_row = 0; F_row < basis_size; ++F_row) {
	  // Check whether this row has already been sampled.
	  std::set<int>::const_iterator found = proc_sampled_f_row[myid].find(F_row);
	  if (found == proc_sampled_f_row[myid].end()) // not found
	    {
	      double tmp = 0.0;
	      for (int F_col = 0; F_col < i; ++F_col) {
		tmp += f_basis->item(F_row, F_col)*c[F_col];
	      }
	      const double r_val = fabs(f_basis->item(F_row, i) - tmp);
	    
	      if (r_val > f_bv_max_local.row_val) {
		f_bv_max_local.row_val = r_val;
		f_bv_max_local.row = F_row;
	      }
	    }
	}
	MPI_Allreduce(&f_bv_max_local, &f_bv_max_global, 1,
		      MaxRowType, RowInfoOp, MPI_COMM_WORLD);

	// Now get the next sampled row of the basis of f.
	if (f_bv_max_global.proc == myid) {
	  for (int j = 0; j < num_basis_vectors; ++j) {
	    c[j] = f_basis->item(f_bv_max_global.row, j);
	  }
	}
	MPI_Bcast(c, num_basis_vectors, MPI_DOUBLE,
		  f_bv_max_global.proc, MPI_COMM_WORLD);
	// Now add the ith sampled row of the basis of the RHS to tmp_fs.
	for (int j = 0; j < num_basis_vectors; ++j) {
	  tmp_fs.item(ns+k, j) = c[j];
	}
	proc_sampled_f_row[f_bv_max_global.proc].insert(f_bv_max_global.row);
	proc_f_row_to_tmp_fs_row[f_bv_max_global.proc][f_bv_max_global.row] = ns+k;
      }

    ns += nsi;
  }

  CAROM_ASSERT(num_samples == ns);
  
  // Fill f_sampled_row, and f_sampled_rows_per_proc.  Unscramble tmp_fs into
  // f_basis_sampled_inv.
  int idx = 0;
  for (int i = 0; i < num_procs; ++i) {
    std::set<int>& this_proc_sampled_f_row = proc_sampled_f_row[i];
    std::map<int, int>& this_proc_f_row_to_tmp_fs_row =
      proc_f_row_to_tmp_fs_row[i];
    f_sampled_rows_per_proc[i] = this_proc_sampled_f_row.size();
    for (std::set<int>::iterator j = this_proc_sampled_f_row.begin();
	 j != this_proc_sampled_f_row.end(); ++j) {
      int this_f_row = *j;
      f_sampled_row[idx] = this_f_row;
      int tmp_fs_row = this_proc_f_row_to_tmp_fs_row[this_f_row];
      for (int col = 0; col < num_f_basis_cols; ++col) {
	f_basis_sampled_inv.item(idx, col) = tmp_fs.item(tmp_fs_row, col);
      }
      ++idx;
    }
  }
  
  CAROM_ASSERT(num_samples == idx);

  // Now invert f_basis_sampled_inv, storing its transpose.
  f_basis_sampled_inv.pseudoinverse();

  // Free the MPI_Datatype and MPI_Op.
  MPI_Type_free(&MaxRowType);
  MPI_Op_free(&RowInfoOp);

  delete [] c;
}

}
