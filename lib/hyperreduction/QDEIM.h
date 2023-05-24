/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Interface to the DEIM algorithm to determine the rows of the
// rhs to be sampled for the interpolation of the rhs.

#ifndef included_QDEIM_h
#define included_QDEIM_h

#include <vector>

namespace CAROM {

class Matrix;

/**
 * @brief Computes the QDEIM algorithm on the given basis.
 *
 * Implemented from Zlatko Drmac, Serkan Gugercin, "A new selection operator for
 * the discrete empirical interpolation method -- Improved a priori error bound
 * and extensions", SIAM Journal on Scientific Computing, Volume 38, Issue 2,
 * pages A631-A648, 2016.
 *
 * Oversampling uses GappyPOD+E from Peherstorfer,
 * Drmac, Gugercin,  "Stability of discrete empirical interpolation and gappy
 * proper orthogonal decomposition with randomized and deterministic sampling
 * points", preprint May 20, 2020.
 *
 * @param[in] f_basis The basis vectors for the RHS.
 * @param[in] num_f_basis_vectors_used The number of basis vectors in f_basis
 *                                     to use in the algorithm.
 * @param[out] f_sampled_row The local row ids of each sampled row.  This will
 *                           contain the sampled rows from all processors.  Use
 *                           f_sampled_rows_per_proc to index into this array
 *                           for the sampled rows from a specific processor.
 * @param[out] f_sampled_rows_per_proc The number of sampled rows for each
 *                                     processor.
 * @param[out] f_basis_sampled_inv The inverse of the sampled basis of the RHS.
 * @param[in] myid The rank of this process.
 * @param[in] num_procs The total number of processes.
 * @param[in] num_samples_req The minimum number of samples required.
 */
void
QDEIM(const Matrix* f_basis,
      int num_f_basis_vectors_used,
      std::vector<int>& f_sampled_row,
      std::vector<int>& f_sampled_rows_per_proc,
      Matrix& f_basis_sampled_inv,
      const int myid,
      const int num_procs,
      const int num_samples_req);
}

#endif
