/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Interface to the DEIM algorithm to determine the rows of the
// rhs to be sampled for the interpolation of the rhs.

#ifndef included_DEIM_h
#define included_DEIM_h

namespace CAROM {

class Matrix;

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
void RowInfoMax(RowInfo* a, RowInfo* b, int* len, MPI_Datatype* type);

/**
 * @brief
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
 * @param[int] myid The rank of this process.
 * @param[int] num_procs The total number of processes.
 */
void
DEIM(const Matrix* f_basis,
     int num_f_basis_vectors_used,
     int* f_sampled_row,
     int* f_sampled_rows_per_proc,
     Matrix& f_basis_sampled_inv,
     int myid,
     int num_procs);

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
GNAT(const Matrix* f_basis,
     const int num_f_basis_vectors_used,
     int* f_sampled_row,
     int* f_sampled_rows_per_proc,
     Matrix& f_basis_sampled_inv,
     const int myid,
     const int num_procs,
     const int num_samples_req = -1);
}

#endif
