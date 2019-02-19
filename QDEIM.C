/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC.
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
