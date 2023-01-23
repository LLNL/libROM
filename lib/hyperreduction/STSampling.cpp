/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

#include "linalg/Matrix.h"
#include "Utilities.h"
#include "mpi.h"
#include <cmath>
#include <map>
#include <set>

#include "STSampling.h"

using namespace std;

namespace CAROM {

// This function implements Algorithm 2 of Choi and Carlberg 2019.
// All spatial indices are used here.
// Spatial basis s_basis is distributed, of size (number of spatial DOFs) x (number of basis vectors)
// Temporal basis t_basis is stored in its entirety on each process, of size (number of temporal DOFs) x (number of basis vectors)
// In general, there may be multiple temporal basis vectors associated with each spatial basis vector, but here we assume for now
// that there is only one, i.e. column j of s_basis corresponds to column j of t_basis.
void SampleTemporalIndices(const Matrix* s_basis,
                           const Matrix* t_basis,
                           const int num_f_basis_vectors_used,
                           int* t_samples,
                           const int myid,
                           const int num_procs,
                           const int num_samples_req,
                           const bool excludeFinalTime)
{
    CAROM_VERIFY(t_basis->distributed());

    // Get the number of basis vectors and the size of each basis vector.
    CAROM_VERIFY(0 < num_f_basis_vectors_used
                 && num_f_basis_vectors_used <= s_basis->numColumns()
                 && num_f_basis_vectors_used <= t_basis->numColumns());
    CAROM_VERIFY(num_samples_req > 0);
    const int num_basis_vectors =
        std::min(num_f_basis_vectors_used, s_basis->numColumns());
    const int num_samples = num_samples_req;
    const int numExcluded = excludeFinalTime ? 1 : 0;
    CAROM_VERIFY(num_samples <= t_basis->numRows() - numExcluded);
    const int s_size = s_basis->numRows();
    const int t_size = t_basis->numRows();

    const int ns_mod_nr = num_samples % num_basis_vectors;
    int ns = 0;

    // The small matrix inverted by the algorithm.  We'll allocate the largest
    // matrix we'll need and set its size at each step in the algorithm.
    Matrix M(num_basis_vectors, num_basis_vectors,
             false);  // TODO: is this big enough to avoid reallocations?

    std::set<int> samples;  // Temporal samples, identical on all processes

    std::vector<double> error(s_size * t_size);
    std::vector<double> MZphi(num_basis_vectors);

    // Initialize error as the first space-time basis vector
    for (int s=0; s<s_size; ++s)
    {
        for (int t=0; t<t_size; ++t)
            error[t + (s*t_size)] = s_basis->item(s, 0) * t_basis->item(t, 0);
    }

    for (int i=0; i<num_basis_vectors; ++i)
    {
        CAROM_VERIFY(samples.size() == ns);
        if (i > 0)
        {
            // Set error vector to (I - [\phi_0, ..., \phi_{i-1}] (Z [\phi_0, ..., \phi_{i-1}])^+ Z) \phi_i
            // where \phi_j is the j-th space-time basis vector (tensor product of columns j of s_basis and t_basis).

            // First, set M to Z [\phi_0, ..., \phi_{i-1}], where Z selects all spatial indices and the temporal indices in `samples`
            M.setSize(ns * s_size, i);

            int ti = 0;
            for (auto t : samples) // loop over ns values
            {
                for (int s = 0; s < s_size; ++s)
                {
                    const int row = ti + (s*ns);
                    for (int col = 0; col < i; ++col)
                    {
                        M.item(row, col) = s_basis->item(s, col) * t_basis->item(t, col);
                    }
                }

                ti++;
            }

            // Compute the pseudo-inverse of M, storing its transpose.
            M.transposePseudoinverse();

            // TODO: in parallel, it seems this matrix should be distributed, and the pseudo-inverse needs to be computed in parallel.
            // In DEIM and GNAT, we gather the non-local rows that have been sampled, but here we use all spatial rows, so this will be
            // more expensive. Fortunately, the matrix to be inverted in the pseudo-inverse computation will still be small.
            // Hopefully M just needs to be declared as `distributed`.

            // Multiply Z\phi_i by the transpose of M, which is (Z [\phi_0, ..., \phi_{i-1}])^+.
            // The result is a vector of length i, stored in MZphi.
            for (int j = 0; j < i; ++j) MZphi[j] = 0.0;

            ti = 0;
            for (auto t : samples) // loop over ns values
            {
                for (int s = 0; s < s_size; ++s)
                {
                    const int row = ti + (s*ns);
                    const double phi_i_s_t = s_basis->item(s, i) * t_basis->item(t,
                                             i);  // \phi_i(s,t)

                    for (int col = 0; col < i; ++col)
                    {
                        MZphi[col] += M.item(row, col) * phi_i_s_t;
                    }
                }

                ti++;
            }

            // Initialize error as space-time basis vector i
            for (int s=0; s<s_size; ++s)
            {
                for (int t=0; t<t_size; ++t)
                    error[t + (s*t_size)] = s_basis->item(s, i) * t_basis->item(t, i);
            }

            // Multiply MZphi by [\phi_0, ..., \phi_{i-1}], subtracting the result from error
            for (int s=0; s<s_size; ++s)
            {
                for (int t=0; t<t_size; ++t)
                {
                    for (int j=0; j<i; ++j)
                        error[t + (s*t_size)] -= MZphi[j] * s_basis->item(s, j) * t_basis->item(t, j);
                }
            }

            // Now error is set.
        }

        const int nsi = i < ns_mod_nr ? (num_samples / num_basis_vectors) + 1 :
                        num_samples / num_basis_vectors;

        for (int j=0; j<nsi; ++j)
        {
            // Find the temporal index not sampled yet that has the maximum l2 spatial error norm

            double maxNorm = -1.0;
            int tmax = 0;
            for (int t=0; t<t_size - numExcluded; ++t)
            {
                std::set<int>::const_iterator found = samples.find(t);
                if (found != samples.end())
                    continue;

                double norm2 = 0.0;
                for (int s=0; s<s_size; ++s)
                {
                    norm2 += error[t + (s*t_size)] * error[t + (s*t_size)];
                }

                // TODO: this Allreduce is expensive for every t. Can this be improved?
                double globalNorm2 = 0.0;
                MPI_Allreduce(&norm2, &globalNorm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                if (globalNorm2 > maxNorm)
                {
                    maxNorm = globalNorm2;
                    tmax = t;
                }
            }

            samples.insert(tmax);
        }

        ns += nsi;
    }

    CAROM_VERIFY(num_samples == ns);

    int ti = 0;
    for (auto t : samples)
    {
        t_samples[ti] = t;
        ti++;
    }
}

// This function implements Algorithm 3 of Choi and Carlberg 2019.
// Temporal samples are input.
// Spatial basis s_basis is distributed, of size (number of spatial DOFs) x (number of basis vectors)
// Temporal basis t_basis is stored in its entirety on each process, of size (number of temporal DOFs) x (number of basis vectors)
// In general, there may be multiple temporal basis vectors associated with each spatial basis vector, but here we assume for now
// that there is only one, i.e. column j of s_basis corresponds to column j of t_basis.
void SampleSpatialIndices(const Matrix* s_basis,
                          const Matrix* t_basis,
                          const int num_f_basis_vectors_used,
                          const int num_t_samples,
                          int* t_samples,
                          int* s_sampled_row,
                          int* s_sampled_rows_per_proc,
                          Matrix& f_basis_sampled,
                          const int myid,
                          const int num_procs,
                          const int num_samples_req)
{
    // This algorithm determines the rows of f that should be sampled, the
    // processor that owns each sampled row, and fills f_basis_sampled with
    // the sampled rows of the basis of the RHS.

    // Create an MPI_Datatype for the RowInfo struct.
    MPI_Datatype MaxRowType, oldtypes[2];
    int blockcounts[2];
    MPI_Aint offsets[2], extent, lower_bound;
    MPI_Status stat;
    offsets[0] = 0;
    oldtypes[0] = MPI_DOUBLE;
    blockcounts[0] = 1;

    /*
        MPI_Type_extent is deprecated in MPI-2 and removed in
        MPI-3. MPI_Type_extent has been replaced by MPI_Type_get_extent.
    */

#if (MPI_VERSION == 1)
    MPI_Type_extent(MPI_DOUBLE, &extent);
#else
    MPI_Type_get_extent(MPI_DOUBLE, &lower_bound, &extent);
#endif

    offsets[1] = extent;
    oldtypes[1] = MPI_INT;
    blockcounts[1] = 2;

    /*
      MPI_Type_struct is deprecated in MPI-2 and removed in
      MPI-3. MPI_Type_struct has been replaced by MPI_Type_create_struct.
     */

#if (MPI_VERSION == 1)
    MPI_Type_struct(2, blockcounts, offsets, oldtypes, &MaxRowType);
#else
    MPI_Type_create_struct(2, blockcounts, offsets, oldtypes, &MaxRowType);
#endif

    MPI_Type_commit(&MaxRowType);

    // Create an MPI_Op for the RowInfoMax function.
    MPI_Op RowInfoOp;
    MPI_Op_create((MPI_User_function*)RowInfoMax, true, &RowInfoOp);

    //CAROM_VERIFY(!t_basis->distributed());

    // Get the number of basis vectors and the size of each basis vector.
    CAROM_VERIFY(0 < num_f_basis_vectors_used
                 && num_f_basis_vectors_used <= s_basis->numColumns()
                 && num_f_basis_vectors_used <= t_basis->numColumns());
    CAROM_VERIFY(num_samples_req >
                 0);  // TODO: just replace num_samples with num_samples_req
    const int num_basis_vectors =
        std::min(num_f_basis_vectors_used, s_basis->numColumns());
    //const int num_samples = num_samples_req > 0 ? num_samples_req : num_basis_vectors;
    const int num_samples = num_samples_req;
    CAROM_VERIFY(num_basis_vectors <= num_samples
                 && num_samples <= s_basis->numRows());
    CAROM_VERIFY(num_samples == f_basis_sampled.numRows()
                 && num_basis_vectors == f_basis_sampled.numColumns());
    CAROM_VERIFY(!f_basis_sampled.distributed());
    const int s_size = s_basis->numRows();
    const int t_size = t_basis->numRows();

    const int ns_mod_nr = num_samples % num_basis_vectors;
    int ns = 0;

    // The small matrix inverted by the algorithm.  We'll allocate the largest
    // matrix we'll need and set its size at each step in the algorithm.
    Matrix M(num_basis_vectors, num_basis_vectors,
             false);  // TODO: is this big enough to avoid reallocations?

    // TODO: change these variable names
    std::vector<std::set<int> > proc_sampled_f_row(num_procs);
    std::vector<std::map<int, int> > proc_f_row_to_tmp_fs_row(num_procs);
    int num_f_basis_cols =
        f_basis_sampled.numColumns(); // TODO: just use num_basis_vectors
    Matrix tmp_fs(f_basis_sampled.numRows(),
                  num_f_basis_cols,
                  f_basis_sampled.distributed()); // TODO: should this be distributed?

    // TODO: change these variable names
    RowInfo f_bv_max_local, f_bv_max_global;

    std::vector<double> error(s_size * t_size);
    std::vector<double> MZphi(num_basis_vectors);
    std::vector<double> sampled_row(num_basis_vectors);

    // Initialize error as the first space-time basis vector
    for (int s=0; s<s_size; ++s)
    {
        for (int t=0; t<t_size; ++t)
            error[t + (s*t_size)] = s_basis->item(s, 0) * t_basis->item(t, 0);
    }

    for (int i=0; i<num_basis_vectors; ++i)
    {
        if (i > 0)
        {
            // Set error vector to (I - [\phi_0, ..., \phi_{i-1}] (Z [\phi_0, ..., \phi_{i-1}])^+ Z) \phi_i
            // where \phi_j is the j-th space-time basis vector (tensor product of columns j of s_basis and t_basis).

            // First, set M to Z [\phi_0, ..., \phi_{i-1}], where Z selects spatial indices in `samples` and the temporal indices in `t_samples`

            /* // Distributed version: should this be used?
            M.setSize(proc_sampled_f_row[myid].size() * num_t_samples, i);

            int si = 0;
            for (auto s : proc_sampled_f_row[myid]) // loop over local spatial sample indices
              {
                for (int ti = 0; ti < num_t_samples; ++ti)
            {
              const int t = t_samples[ti];
              const int row = ti + (si*num_t_samples);
              for (int col = 0; col < i; ++col)
                {
                  M.item(row, col) = s_basis->item(s, col) * t_basis->item(t, col);
                }
            }

                si++;
              }
            */

            // Set M to be the global sampled matrix and invert on each process.
            M.setSize(ns * num_t_samples, i);
            for (int k = 0; k < i; ++k)
                for (int j = 0; j < ns; ++j)
                {
                    for (int ti = 0; ti < num_t_samples; ++ti)
                    {
                        const int t = t_samples[ti];
                        const int row = ti + (j*num_t_samples);

                        M.item(row, k) = tmp_fs.item(j, k) * t_basis->item(t, k);
                    }
                }

            // Compute the pseudo-inverse of M, storing its transpose.
            M.transposePseudoinverse();

            // Multiply Z\phi_i by the transpose of M, which is (Z [\phi_0, ..., \phi_{i-1}])^+.
            // The result is a vector of length i, stored in MZphi.
            for (int j = 0; j < i; ++j) MZphi[j] = 0.0;

            /*
            // Distributed version: should this be used?

            // TODO: in parallel, it seems this matrix should be distributed, and the pseudo-inverse needs to be computed in parallel.
            // In DEIM and GNAT, we gather the non-local spatial rows that have been sampled. This is so far just a serial implementation.

            si = 0;
            for (auto s : proc_sampled_f_row[myid]) // loop over local spatial sample indices
              {
                for (int ti = 0; ti < num_t_samples; ++ti)
            {
              const int t = t_samples[ti];
              const int row = ti + (si*num_t_samples);
              const double phi_i_s_t = s_basis->item(s, i) * t_basis->item(t, i);  // \phi_i(s,t)
              for (int col = 0; col < i; ++col)
                {
                  MZphi[col] += M.item(row, col) * phi_i_s_t;
                }
            }

                si++;
              }

            CAROM_VERIFY(si * num_t_samples == M.numRows());
            */

            for (int j = 0; j < ns; ++j)
            {
                for (int ti = 0; ti < num_t_samples; ++ti)
                {
                    const int t = t_samples[ti];
                    const int row = ti + (j*num_t_samples);

                    const double phi_i_s_t = s_basis->item(j, i) * t_basis->item(t,
                                             i);  // \phi_i(s,t)

                    for (int k = 0; k < i; ++k)
                        MZphi[k] += M.item(row, k) * phi_i_s_t;
                }
            }

            // Initialize error as space-time basis vector i
            for (int s=0; s<s_size; ++s)
            {
                for (int t=0; t<t_size; ++t)
                    error[t + (s*t_size)] = s_basis->item(s, i) * t_basis->item(t, i);
            }

            // Multiply MZphi by [\phi_0, ..., \phi_{i-1}], subtracting the result from error
            for (int s=0; s<s_size; ++s)
            {
                for (int t=0; t<t_size; ++t)
                {
                    for (int j=0; j<i; ++j)
                        error[t + (s*t_size)] -= MZphi[j] * s_basis->item(s, j) * t_basis->item(t, j);
                }
            }

            // Now error is set.
        }

        const int nsi = i < ns_mod_nr ? (num_samples / num_basis_vectors) + 1 :
                        num_samples / num_basis_vectors;

        for (int j=0; j<nsi; ++j)
        {
            // Find the spatial index not sampled yet that has the maximum l2 temporal error norm.
            // First find the process-local maximum.

            double maxNorm = -1.0;
            int smax = 0;
            for (int s=0; s<s_size; ++s)
            {
                std::set<int>::const_iterator found = proc_sampled_f_row[myid].find(s);
                if (found != proc_sampled_f_row[myid].end())
                    continue;

                double norm2 = 0.0;
                for (int t=0; t<t_size; ++t)
                {
                    norm2 += error[t + (s*t_size)] * error[t + (s*t_size)];
                }

                if (norm2 > maxNorm)
                {
                    maxNorm = norm2;
                    smax = s;
                }
            }

            // The local maximum with respect spatial indices is (maxNorm, smax).
            // Now find the global maximum.
            //RowInfo f_bv_max_local, f_bv_max_global;
            f_bv_max_local.row_val = maxNorm;
            f_bv_max_local.row = smax;
            f_bv_max_local.proc = myid;

            MPI_Allreduce(&f_bv_max_local, &f_bv_max_global, 1,
                          MaxRowType, RowInfoOp, MPI_COMM_WORLD);

            proc_sampled_f_row[f_bv_max_global.proc].insert(f_bv_max_global.row);
            proc_f_row_to_tmp_fs_row[f_bv_max_global.proc][f_bv_max_global.row] = ns + j;

            // Now get the sampled row of the spatial basis.
            if (f_bv_max_global.proc == myid) {
                for (int k = 0; k < num_basis_vectors; ++k) {
                    sampled_row[k] = s_basis->item(f_bv_max_global.row, k);
                }
            }
            MPI_Bcast(sampled_row.data(), num_basis_vectors, MPI_DOUBLE,
                      f_bv_max_global.proc, MPI_COMM_WORLD);
            // Now add the sampled row of the basis to tmp_fs.
            for (int k = 0; k < num_basis_vectors; ++k) {
                tmp_fs.item(ns+j, k) = sampled_row[k];
            }
        }

        ns += nsi;
    }

    CAROM_VERIFY(num_samples == ns);

    // Fill f_sampled_row, and f_sampled_rows_per_proc.  Unscramble tmp_fs into
    // f_basis_sampled.
    int idx = 0;
    for (int i = 0; i < num_procs; ++i) {
        std::set<int>& this_proc_sampled_f_row = proc_sampled_f_row[i];
        std::map<int, int>& this_proc_f_row_to_tmp_fs_row =
            proc_f_row_to_tmp_fs_row[i];
        s_sampled_rows_per_proc[i] = this_proc_sampled_f_row.size();
        for (std::set<int>::iterator j = this_proc_sampled_f_row.begin();
                j != this_proc_sampled_f_row.end(); ++j) {
            int this_f_row = *j;
            s_sampled_row[idx] = this_f_row;
            int tmp_fs_row = this_proc_f_row_to_tmp_fs_row[this_f_row];
            for (int col = 0; col < num_f_basis_cols; ++col) {
                f_basis_sampled.item(idx, col) = tmp_fs.item(tmp_fs_row, col);
            }
            ++idx;
        }
    }

    CAROM_VERIFY(num_samples == idx);

    // Free the MPI_Datatype and MPI_Op.
    MPI_Type_free(&MaxRowType);
    MPI_Op_free(&RowInfoOp);
}

void SpaceTimeSampling(const Matrix* s_basis,
                       const Matrix* t_basis,
                       const int num_f_basis_vectors_used,
                       std::vector<int>& t_samples,
                       int* f_sampled_row,
                       int* f_sampled_rows_per_proc,
                       Matrix& s_basis_sampled,
                       const int myid,
                       const int num_procs,
                       const int num_t_samples_req,
                       const int num_s_samples_req,
                       const bool excludeFinalTime)
{
    // There are multiple possible algorithms for space-time sampling. For now, we just implement
    // one algorithm, but in the future there could be options for other algorithms.
    // This algorithm is sequential greedy sampling of temporal then spatial indices.

    // TODO: for now, we assume one temporal basis vector for each spatial basis vector. This should be generalized.
    CAROM_VERIFY(s_basis->numColumns() == num_f_basis_vectors_used
                 && t_basis->numColumns() == num_f_basis_vectors_used);

    // First, sample temporal indices.
    CAROM_VERIFY(t_samples.size() == num_t_samples_req);

    SampleTemporalIndices(s_basis, t_basis, num_f_basis_vectors_used,
                          t_samples.data(),
                          myid, num_procs, num_t_samples_req, excludeFinalTime);

    // Second, sample spatial indices.
    //Matrix s_basis_sampled(num_s_samples_req, num_f_basis_vectors_used, false);
    CAROM_VERIFY(s_basis_sampled.numRows() == num_s_samples_req
                 && s_basis_sampled.numColumns() == num_f_basis_vectors_used);

    //std::vector<int> s_samples(num_s_samples_req);
    //std::vector<int> s_samples_per_proc(num_procs);
    SampleSpatialIndices(s_basis, t_basis, num_f_basis_vectors_used,
                         num_t_samples_req,
                         //t_samples.data(), s_samples.data(), s_samples_per_proc.data(),
                         t_samples.data(), f_sampled_row, f_sampled_rows_per_proc,
                         s_basis_sampled, myid, num_procs, num_s_samples_req);
}

void GetSampledSpaceTimeBasis(std::vector<int> const& t_samples,
                              const Matrix* t_basis,
                              Matrix const& s_basis_sampled,
                              Matrix& f_basis_sampled_inv)
{
    const int num_s_samples = s_basis_sampled.numRows();
    const int num_t_samples = t_samples.size();

    // Set sampled space-time basis matrix and take its pseudo-inverse, in f_basis_sampled_inv.

    CAROM_VERIFY(f_basis_sampled_inv.numRows() == num_t_samples * num_s_samples);

    for (int si=0; si<num_s_samples; ++si)
    {
        for (int ti=0; ti<num_t_samples; ++ti)
        {
            const int row = ti + (si*num_t_samples);
            const int t = t_samples[ti];
            for (int j=0; j<f_basis_sampled_inv.numColumns(); ++j)
                f_basis_sampled_inv.item(row, j) = s_basis_sampled.item(si,
                                                   j) * t_basis->item(t, j);
        }
    }

    // Compute the pseudo-inverse of f_basis_sampled_inv, storing its transpose.
    f_basis_sampled_inv.transposePseudoinverse();
}

}
