/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Interface to the DEIM algorithm to determine the rows of the
// rhs to be sampled for the interpolation of the rhs.

#include "Matrix.h"
#include "mpi.h"
#include <cmath>

namespace CAROM {

  /* Implement QDEIM algorithm from Zlatko Drmac, Serkan Gugercin, "A
     new selection operator for the discrete empirical interpolation
     method -- Improved a priori error bound and extensions", SIAM
     Journal on Scientific Computing, Volume 38, Issue 2, pages
     A631-A648, 2016.
  */
void
QDEIM(const Matrix* f_basis,
      int num_f_basis_vectors_used,
      int* f_sampled_row,
      int* f_sampled_rows_per_proc,
      Matrix& f_basis_sampled_inv,
      const int myid,
      const int num_procs,
      const int num_samples_req)
{
  // This algorithm determines the rows of f that should be sampled, the
  // processor that owns each sampled row, and fills f_basis_sampled_inv with the
  // sampled rows of the basis of the RHS.

  CAROM_VERIFY(num_f_basis_vectors_used == f_basis->numColumns());  // The QR implementation uses the entire matrix.
  CAROM_VERIFY(f_basis->numColumns() <= num_samples_req && num_samples_req <= f_basis->numRows());
  
  int *f_sampled_row_owner = new int[num_samples_req];

  // QDEIM computes selection/interpolation indices by taking a
  // column-pivoted QR-decomposition of the transpose of its input
  // matrix.
  f_basis->qrcp_pivots_transpose(f_sampled_row,
				 f_sampled_row_owner,
				 num_samples_req);

  if (f_basis->distributed())
    {
      // Gather the sampled rows to root process in f_basis_sampled_inv

      // On root, f_sampled_row* contains all global pivots.
      // On non-root processes, they contain the local pivots.
      
      int * n = (myid == 0) ? new int[num_procs] : NULL;
      int * disp = (myid == 0) ? new int[num_procs] : NULL;
      int * all_sampled_rows = (myid == 0) ? new int[num_samples_req] : NULL;
      if (myid == 0)
	{
	  for (int r=0; r<num_procs; ++r)
	    n[r] = 0;

	  for (int i=0; i<num_samples_req; ++i)
	    n[f_sampled_row_owner[i]]++;

	  disp[0] = 0;
	  for (int r=1; r<num_procs; ++r)
	    disp[r] = disp[r-1] + n[r-1];

	  CAROM_VERIFY(disp[num_procs-1] + n[num_procs-1] == num_samples_req);

	  for (int r=0; r<num_procs; ++r)
	    {
	      f_sampled_rows_per_proc[r] = n[r];
	      n[r] = 0;
	    }
	  
	  for (int i=0; i<num_samples_req; ++i)
	    {
	      const int owner = f_sampled_row_owner[i];
	      all_sampled_rows[disp[owner] + n[owner]] = f_sampled_row[i];
	      n[owner]++;
	    }
	}

      MPI_Bcast(f_sampled_rows_per_proc, num_procs, MPI_INT, 0, MPI_COMM_WORLD);
      
      // TODO: eliminate n and the following MPI_Scatter
      
      int count = 0;
      MPI_Scatter(n, 1, MPI_INT, &count, 1, MPI_INT, 0, MPI_COMM_WORLD);

      int * my_sampled_rows = (count > 0) ? new int[count] : NULL;
      double * my_sampled_row_data = (count > 0) ? new double[count*num_f_basis_vectors_used] : NULL;

      MPI_Scatterv(all_sampled_rows, n, disp, MPI_INT, my_sampled_rows, count, MPI_INT, 0, MPI_COMM_WORLD);

      int *row_offset = new int[num_procs];
      row_offset[myid] = f_basis->numRows();

      CAROM_VERIFY(MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, row_offset, 1,
				 MPI_INT, MPI_COMM_WORLD) == MPI_SUCCESS);

      int os = 0;
      for (int i=0; i<num_procs; ++i)
	{
	  os += row_offset[i];
	  row_offset[i] = os - row_offset[i];
	}
      
      for (int i=0; i<count; ++i)
	{
	  CAROM_VERIFY(my_sampled_rows[i] >= row_offset[myid] && my_sampled_rows[i] < row_offset[myid] + f_basis->numRows());
	  const int row = my_sampled_rows[i] - row_offset[myid];
	  os = i*num_f_basis_vectors_used;
	  for (int j=0; j<num_f_basis_vectors_used; ++j)
	    my_sampled_row_data[os + j] = f_basis->item(row, j);
	}

      if (myid == 0)
	{
	  for (int r=0; r<num_procs; ++r)
	    {
	      n[r] *= num_f_basis_vectors_used;
	      disp[r] *= num_f_basis_vectors_used;
	    }
	}
      
      MPI_Gatherv(my_sampled_row_data, count*num_f_basis_vectors_used, MPI_DOUBLE, f_basis_sampled_inv.getData(), n, disp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      
      // Subtract row_offset to convert f_sampled_row from global to local indices
      os = row_offset[myid];
      for (int i=0; i<num_f_basis_vectors_used; ++i)
	{
	  f_sampled_row[i] -= os;
	}

      delete [] n;
      delete [] disp;
      delete [] my_sampled_rows;
      delete [] all_sampled_rows;
      delete [] row_offset;
    }
  else
    {
      // With the known interpolation (sample) indices, copy over the
      // rows of the sampled basis
      for (int i = 0; i < num_samples_req; i++) {
	for (int j = 0; j < num_f_basis_vectors_used; j++) {
	  f_basis_sampled_inv.item(i, j) = f_basis->item(f_sampled_row[i], j);
	}
      }

      f_sampled_rows_per_proc[0] = num_f_basis_vectors_used;
    }

  delete [] f_sampled_row_owner;
  
  // Now invert f_basis_sampled_inv, storing its transpose.
  f_basis_sampled_inv.transposePseudoinverse();
} // end void QDEIM

} // end namespace CAROM
