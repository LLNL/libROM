/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Interface to the space-time sampling index selection algorithm.

#ifndef included_STSAMPLING_h
#define included_STSAMPLING_h

namespace CAROM {

  // TODO: update this documentation
/**
 * @brief
 *
 * @param[in] f_basis The basis vectors for the RHS.
 * @param[in] num_samples The number of samples to compute.
 * @param[in] num_f_basis_vectors_used The number of basis vectors in f_basis
 *                                     to use in the algorithm.
 * @param[out] f_sampled_row The local row ids of each sampled row.  This will
 *                           contain the sampled rows from all processors.  Use
 *                           f_sampled_rows_per_proc to index into this array
 *                           for the sampled rows from a specific processor.
 * @param[out] f_sampled_rows_per_proc The number of sampled rows for each
 *                                     processor.
 * @param[out] f_basis_sampled_inv The inverse of the sampled basis of the RHS.
 * @param[int] myid The rank of this process.
 * @param[int] num_procs The total number of processes.
 */
void
SpaceTimeSampling(const Matrix* s_basis,
		  const Matrix* t_basis,
		  const int num_f_basis_vectors_used,
		  std::vector<int>& t_samples,
		  int* f_sampled_row,
		  int* f_sampled_rows_per_proc,
		  Matrix& s_basis_sampled,
		  const int myid,
		  const int num_procs,
		  const int num_t_samples_req = -1,  // TODO: remove defaults?
		  const int num_s_samples_req = -1,
		  const bool excludeFinalTime = false);

 void GetSampledSpaceTimeBasis(std::vector<int> const& t_samples,
			       const Matrix* t_basis,
			       Matrix const& s_basis_sampled,
			       Matrix& f_basis_sampled_inv);

}

#endif
