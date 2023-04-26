/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Interface to the GNAT algorithm to determine the rows of the
// rhs to be sampled for the interpolation of the rhs.

#include "linalg/Matrix.h"
#include "Utilities.h"
#include "mpi.h"
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include "GNAT.h"

namespace CAROM {

void GNAT(const Matrix* f_basis,
          const int num_f_basis_vectors_used,
          std::vector<int>& f_sampled_row,
          std::vector<int>& f_sampled_rows_per_proc,
          Matrix& f_basis_sampled_inv,
          const int myid,
          const int num_procs,
          const int num_samples_req,
          std::vector<int> *init_samples)
{
    CAROM_VERIFY(num_procs == f_sampled_rows_per_proc.size());

    // This algorithm determines the rows of f that should be sampled, the
    // processor that owns each sampled row, and fills f_basis_sampled_inv with
    // the inverse of the sampled rows of the basis of the RHS.

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

    // Get the number of basis vectors and the size of each basis vector.
    CAROM_VERIFY(0 < num_f_basis_vectors_used
                 && num_f_basis_vectors_used <= f_basis->numColumns());
    const int num_basis_vectors =
        std::min(num_f_basis_vectors_used, f_basis->numColumns());
    const int num_samples = num_samples_req > 0 ? num_samples_req :
                            num_basis_vectors;
    CAROM_VERIFY(num_basis_vectors <= num_samples
                 && num_samples <= f_basis->numDistributedRows());
    CAROM_VERIFY(num_samples == f_sampled_row.size());
    CAROM_VERIFY(num_samples == f_basis_sampled_inv.numRows()
                 && num_basis_vectors == f_basis_sampled_inv.numColumns());
    CAROM_VERIFY(!f_basis_sampled_inv.distributed());
    const int basis_size = f_basis->numRows();

    const int ns_mod_nr = num_samples % num_basis_vectors;
    int ns = 0;

    // The small matrix inverted by the algorithm.  We'll allocate the largest
    // matrix we'll need and set its size at each step in the algorithm.
    Matrix M(num_samples, std::max(num_basis_vectors-1, 1), false);

    // Scratch space used throughout the algorithm.
    double* c = new double [num_basis_vectors];
    double* sampled_row = new double [num_basis_vectors];

    std::vector<std::set<int> > proc_sampled_f_row(num_procs);
    std::vector<std::map<int, int> > proc_f_row_to_tmp_fs_row(num_procs);
    int num_f_basis_cols = f_basis_sampled_inv.numColumns();
    Matrix tmp_fs(f_basis_sampled_inv.numRows(),
                  num_f_basis_cols,
                  f_basis_sampled_inv.distributed());

    // Gather information about initial samples given as input.
    const int num_init_samples = init_samples ? init_samples->size() : 0;
    int total_num_init_samples = 0;
    MPI_Allreduce(&num_init_samples, &total_num_init_samples, 1,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    int init_sample_offset = 0;
    if (total_num_init_samples > 0)
    {
        CAROM_VERIFY(init_samples);
        std::vector<int> all_num_init_samples(num_procs);
        std::vector<int> all_init_samples(total_num_init_samples);

        MPI_Allgather(&num_init_samples, 1, MPI_INT, all_num_init_samples.data(), 1,
                      MPI_INT, MPI_COMM_WORLD);

        for (int i = 0; i < myid; ++i)
        {
            init_sample_offset += all_num_init_samples[i];
        }
    }

    // Figure out the 1st sampled rows of the RHS.
    RowInfo f_bv_max_local, f_bv_max_global;

    const int ns0 = 0 < ns_mod_nr ? (num_samples / num_basis_vectors) + 1 :
                    num_samples / num_basis_vectors;

    for (int k=0; k<ns0; ++k)
    {
        f_bv_max_local.row_val = -1.0;
        f_bv_max_local.proc = myid;

        if (k < total_num_init_samples)
        {
            // Take sample from init_samples on the corresponding process
            if (k >= init_sample_offset && k - init_sample_offset < num_init_samples)
            {
                f_bv_max_local.row_val = 1.0;  // arbitrary number, ensuring maximum
                f_bv_max_local.row = (*init_samples)[k - init_sample_offset];
            }
        }
        else
        {
            // Compute sample by the greedy algorithm
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
        }

        MPI_Allreduce(&f_bv_max_local, &f_bv_max_global, 1,
                      MaxRowType, RowInfoOp, MPI_COMM_WORLD);

        // Now get the first sampled row of the basis of the RHS.
        if (f_bv_max_global.proc == myid) {
            for (int j = 0; j < num_basis_vectors; ++j) {
                sampled_row[j] = f_basis->item(f_bv_max_global.row, j);
            }
        }
        MPI_Bcast(sampled_row, num_basis_vectors, MPI_DOUBLE,
                  f_bv_max_global.proc, MPI_COMM_WORLD);
        // Now add the first sampled row of the basis of the RHS to tmp_fs.
        for (int j = 0; j < num_basis_vectors; ++j) {
            tmp_fs.item(k, j) = sampled_row[j];
        }
        proc_sampled_f_row[f_bv_max_global.proc].insert(f_bv_max_global.row);
        proc_f_row_to_tmp_fs_row[f_bv_max_global.proc][f_bv_max_global.row] = k;
    }

    ns += ns0;

    // Now repeat the process for the other sampled rows of the basis of the RHS.
    for (int i = 1; i < num_basis_vectors; ++i) {
        const int nsi = i < ns_mod_nr ? (num_samples / num_basis_vectors) + 1 :
                        num_samples / num_basis_vectors;

        // If we currently know about S sampled rows of the basis of the RHS then
        // M contains the first i columns of those S sampled rows.
        M.setSize(ns, i);
        for (int row = 0; row < ns; ++row) {
            for (int col = 0; col < i; ++col) {
                M.item(row, col) = tmp_fs.item(row, col);
            }
        }

        // Compute the pseudo-inverse of M, storing its transpose.
        M.transposePseudoinverse();

        // Now compute c, the inverse of M times the next column of the sampled
        // rows of the basis of the RHS.
        for (int minv_row = 0; minv_row < i; ++minv_row) {
            double tmp = 0.0;
            for (int minv_col = 0; minv_col < ns; ++minv_col) {
                if (ns == i)
                    tmp += M.item(minv_row, minv_col)*tmp_fs.item(minv_col, i);
                else
                    tmp += M.item(minv_col, minv_row)*tmp_fs.item(minv_col,
                            i);  // Transposing M^+, which is stored as its transpose.
            }
            c[minv_row] = tmp;
        }

        for (int k=0; k<nsi; ++k)
        {
            // Now figure out the next sampled row of the basis of f.
            // Compute the first S basis vectors of the RHS times c and find the
            // row of this product have the greatest absolute value. This is the
            // next sampled row of the basis of f.
            f_bv_max_local.row_val = -1.0;
            f_bv_max_local.proc = myid;

            if (ns + k < total_num_init_samples)
            {
                // Take sample from init_samples on the corresponding process
                if (ns + k >= init_sample_offset
                        && ns + k - init_sample_offset < num_init_samples)
                {
                    f_bv_max_local.row_val = 1.0;  // arbitrary number, ensuring maximum
                    f_bv_max_local.row = (*init_samples)[ns + k - init_sample_offset];
                }
            }
            else
            {
                // Compute sample by the greedy algorithm
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
            }

            MPI_Allreduce(&f_bv_max_local, &f_bv_max_global, 1,
                          MaxRowType, RowInfoOp, MPI_COMM_WORLD);

            // Now get the next sampled row of the basis of f.
            if (f_bv_max_global.proc == myid) {
                for (int j = 0; j < num_basis_vectors; ++j) {
                    sampled_row[j] = f_basis->item(f_bv_max_global.row, j);
                }
            }
            MPI_Bcast(sampled_row, num_basis_vectors, MPI_DOUBLE,
                      f_bv_max_global.proc, MPI_COMM_WORLD);
            // Now add the ith sampled row of the basis of the RHS to tmp_fs.
            for (int j = 0; j < num_basis_vectors; ++j) {
                tmp_fs.item(ns+k, j) = sampled_row[j];
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
    f_basis_sampled_inv.transposePseudoinverse();

    // Free the MPI_Datatype and MPI_Op.
    MPI_Type_free(&MaxRowType);
    MPI_Op_free(&RowInfoOp);

    delete [] c;
    delete [] sampled_row;
}

}
