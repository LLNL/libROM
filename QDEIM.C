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
      int myid)
{
  // This algorithm determines the rows of f that should be sampled, the
  // processor that owns each sampled row, and fills f_basis_sampled with the
  // sampled rows of the basis of the RHS.

  // Compute the number of basis vectors used; this number can't
  // exceed the number of columns of the matrix.
  int num_basis_vectors_used = std::min(num_f_basis_vectors_used,
					f_basis->numColumns());

  // QDEIM computes selection/interpolation indices by taking a
  // column-pivoted QR-decomposition of the transpose of its input
  // matrix.
  f_basis->qrcp_pivots_transpose(f_sampled_row,
				 f_sampled_row_owner,
				 num_basis_vectors_used);

  // With the known interpolation (sample) indices, copy over the
  // rows of the sampled basis
  for (int i = 0; i < num_basis_vectors_used; i++) {
    for (int j = 0; j < num_basis_vectors_used; j++) {
      f_basis_sampled.item(i, j) = f_basis->item(f_sampled_row[i], j);
    }
  }

} // end void QDEIM

} // end namespace CAROM
