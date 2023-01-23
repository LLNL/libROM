/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Interface to the S_OPT algorithm to determine the rows of the
// basis to be sampled for the interpolation of the basis.

#ifndef included_S_OPT_h
#define included_S_OPT_h

#include <vector>

namespace CAROM {

class Matrix;

/**
 * @brief Computes the S_OPT algorithm on the given basis.
 *
 * Implemented from Yeonjong Shin and Dongbin Xiu's "Nonadaptive Quasi-Optimal Points
 * Selection for Least Squares Linear Regression", SIAM J. Sci. Comput., 38(1),
 * A385-A411. (26 pages)
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
 * @param[in] init_samples Samples to initialize the S_OPT algorithm.
 * @param[in] qr_factorize Whether to factorize the incoming matrix. If true and
 *                         if the incoming matrix is a basis, the unnecessary
 *                         computation of a QR factorization will be performed.
 */
void
S_OPT(const Matrix* f_basis,
      int num_f_basis_vectors_used,
      std::vector<int>& f_sampled_row,
      std::vector<int>& f_sampled_rows_per_proc,
      Matrix& f_basis_sampled_inv,
      const int myid,
      const int num_procs,
      const int num_samples_req = -1,
      std::vector<int> *init_samples=NULL,
      bool qr_factorize = false);

}

#endif
