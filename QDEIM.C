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
      int* f_sampled_row_owner,
      Matrix& f_basis_sampled,
      const int myid)
{
  // This algorithm determines the rows of f that should be sampled, the
  // processor that owns each sampled row, and fills f_basis_sampled with the
  // sampled rows of the basis of the RHS.

  // Compute the number of basis vectors used; this number can't
  // exceed the number of columns of the matrix.
  const int num_basis_vectors_used = std::min(num_f_basis_vectors_used,
					      f_basis->numColumns());

  // QDEIM computes selection/interpolation indices by taking a
  // column-pivoted QR-decomposition of the transpose of its input
  // matrix.
  f_basis->qrcp_pivots_transpose(f_sampled_row,
				 f_sampled_row_owner,
				 num_basis_vectors_used);

  if (f_basis->distributed())
    {
      // Gather the sampled rows to root process in f_basis_sampled

      // On root, f_sampled_row* contains all global pivots.
      // On non-root processes, they contain the local pivots.

      int rank, num_procs;
      {
	const bool success = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	CAROM_ASSERT(success);
	CAROM_ASSERT(rank == myid);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
      }
      
      int * n = (rank == 0) ? new int[num_procs] : NULL;
      int * disp = (rank == 0) ? new int[num_procs] : NULL;
      int * all_sampled_rows = (rank == 0) ? new int[num_basis_vectors_used] : NULL;
      if (rank == 0)
	{
	  for (int r=0; r<num_procs; ++r)
	    n[r] = 0;

	  for (int i=0; i<num_basis_vectors_used; ++i)
	    n[f_sampled_row_owner[i]]++;

	  disp[0] = 0;
	  for (int r=1; r<num_procs; ++r)
	    disp[r] = disp[r-1] + n[r-1];

	  CAROM_ASSERT(disp[num_procs-1] + n[num_procs-1] == num_basis_vectors_used);

	  for (int r=0; r<num_procs; ++r)
	    n[r] = 0;

	  for (int i=0; i<num_basis_vectors_used; ++i)
	    {
	      const int owner = f_sampled_row_owner[i];
	      all_sampled_rows[disp[owner] + n[owner]] = f_sampled_row[i];
	      n[owner]++;
	    }
	}

      int count = 0;
      MPI_Scatter(n, 1, MPI_INT, &count, 1, MPI_INT, 0, MPI_COMM_WORLD);

      int * my_sampled_rows = (count > 0) ? new int[count] : NULL;
      double * my_sampled_row_data = (count > 0) ? new double[count*num_basis_vectors_used] : NULL;

      MPI_Scatterv(all_sampled_rows, n, disp, MPI_INT, my_sampled_rows, count, MPI_INT, 0, MPI_COMM_WORLD);

      int *row_offset = new int[num_procs];
      row_offset[rank] = f_basis->numRows();

      CAROM_ASSERT(MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, row_offset, 1,
				 MPI_INT, MPI_COMM_WORLD) == MPI_SUCCESS);

      int os = 0;
      for (int i=0; i<rank; ++i)
	{
	  os += row_offset[i];
	  row_offset[i] = os - row_offset[i];
	}
      
      for (int i=0; i<count; ++i)
	{
	  CAROM_ASSERT(my_sampled_rows[i] >= row_offset[rank] && my_sampled_rows[i] < row_offset[rank] + f_basis->numRows());
	  const int row = my_sampled_rows[i] - row_offset[rank];
	  os = i*num_basis_vectors_used;
	  for (int j=0; j<num_basis_vectors_used; ++j)
	    my_sampled_row_data[os + j] = f_basis->item(row, j);
	}

      if (rank == 0)
	{
	  for (int r=0; r<num_procs; ++r)
	    {
	      n[r] *= num_basis_vectors_used;
	      disp[r] *= num_basis_vectors_used;
	    }
	}
      
      MPI_Gatherv(my_sampled_row_data, count*num_basis_vectors_used, MPI_DOUBLE, f_basis_sampled.getData(), n, disp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      
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
      for (int i = 0; i < num_basis_vectors_used; i++) {
	for (int j = 0; j < num_basis_vectors_used; j++) {
	  f_basis_sampled.item(i, j) = f_basis->item(f_sampled_row[i], j);
	}
      }
    }
  
} // end void QDEIM

} // end namespace CAROM
