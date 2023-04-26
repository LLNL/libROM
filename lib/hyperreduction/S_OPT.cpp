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

#include "linalg/Matrix.h"
#include "Utilities.h"
#include "mpi.h"
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <cfloat>

#include "S_OPT.h"

using namespace std;

namespace CAROM {

void
S_OPT(const Matrix* f_basis,
      int num_f_basis_vectors_used,
      std::vector<int>& f_sampled_row,
      std::vector<int>& f_sampled_rows_per_proc,
      Matrix& f_basis_sampled_inv,
      const int myid,
      const int num_procs,
      const int num_samples_req,
      std::vector<int>* init_samples,
      bool qr_factorize)
{
    CAROM_VERIFY(num_procs == f_sampled_rows_per_proc.size());
    // This algorithm determines the rows of f that should be sampled, the
    // processor that owns each sampled row, and fills f_basis_sampled_inv with
    // the inverse of the sampled rows of the basis.

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
    CAROM_VERIFY(num_samples == f_basis_sampled_inv.numRows() &&
                 num_basis_vectors == f_basis_sampled_inv.numColumns());
    CAROM_VERIFY(!f_basis_sampled_inv.distributed());

    int num_rows = f_basis->numRows();

    // If num_basis_vectors is less than the number of columns of the basis,
    // we need to truncate the basis.
    const Matrix* f_basis_truncated = NULL;
    if (num_basis_vectors < f_basis->numColumns())
    {
        f_basis_truncated = f_basis->getFirstNColumns(num_basis_vectors);
    }
    else
    {
        f_basis_truncated = f_basis;
    }

    const Matrix* Vo = NULL;

    // Use the QR factorization of the input matrix, if requested
    if (qr_factorize)
    {
        Vo = f_basis_truncated->qr_factorize();
    }
    else
    {
        Vo = f_basis_truncated;
    }

    int num_samples_obtained = 0;

    // Scratch space used throughout the algorithm.
    std::vector<double> c(num_basis_vectors);

    vector<set<int> > proc_sampled_f_row(num_procs);
    vector<map<int, int> > proc_f_row_to_tmp_fs_row(num_procs);
    int num_f_basis_cols = f_basis_sampled_inv.numColumns();
    Matrix V1(num_samples,
              num_basis_vectors,
              false);

    RowInfo f_bv_max_local, f_bv_max_global;

    // Gather information about initial samples given as input.
    const int num_init_samples = init_samples ? init_samples->size() : 0;
    int total_num_init_samples = 0;
    MPI_Allreduce(&num_init_samples, &total_num_init_samples, 1,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CAROM_VERIFY(num_samples >= total_num_init_samples);

    int init_sample_offset = 0;
    for (int i = 0; i < total_num_init_samples; i++)
    {
        f_bv_max_local.row_val = -DBL_MAX;
        f_bv_max_local.proc = myid;
        if (init_sample_offset < num_init_samples)
        {
            f_bv_max_local.row_val = 1.0;
            f_bv_max_local.row = (*init_samples)[init_sample_offset];
        }
        MPI_Allreduce(&f_bv_max_local, &f_bv_max_global, 1,
                      MaxRowType, RowInfoOp, MPI_COMM_WORLD);
        // Now get the first sampled row of the basis
        if (f_bv_max_global.proc == myid) {
            for (int j = 0; j < num_basis_vectors; ++j) {
                c[j] = Vo->item(f_bv_max_global.row, j);
            }
            init_sample_offset++;
        }
        MPI_Bcast(c.data(), num_basis_vectors, MPI_DOUBLE,
                  f_bv_max_global.proc, MPI_COMM_WORLD);
        // Now add the first sampled row of the basis to tmp_fs.
        for (int j = 0; j < num_basis_vectors; ++j) {
            V1.item(num_samples_obtained, j) = c[j];
        }
        proc_sampled_f_row[f_bv_max_global.proc].insert(f_bv_max_global.row);
        proc_f_row_to_tmp_fs_row[f_bv_max_global.proc][f_bv_max_global.row] =
            num_samples_obtained;
        num_samples_obtained++;
    }

    if (num_samples_obtained == 0)
    {
        f_bv_max_local.row_val = -DBL_MAX;
        f_bv_max_local.proc = myid;
        for (int i = 0; i < num_rows; ++i) {
            double f_bv_val = fabs(Vo->item(i, 0));
            if (f_bv_val > f_bv_max_local.row_val) {
                f_bv_max_local.row_val = f_bv_val;
                f_bv_max_local.row = i;
            }
        }
        MPI_Allreduce(&f_bv_max_local, &f_bv_max_global, 1,
                      MaxRowType, RowInfoOp, MPI_COMM_WORLD);
        // Now get the first sampled row of the basis
        if (f_bv_max_global.proc == myid) {
            for (int j = 0; j < num_basis_vectors; ++j) {
                c[j] = Vo->item(f_bv_max_global.row, j);
            }
        }
        MPI_Bcast(c.data(), num_basis_vectors, MPI_DOUBLE,
                  f_bv_max_global.proc, MPI_COMM_WORLD);
        // Now add the first sampled row of the basis to tmp_fs.
        for (int j = 0; j < num_basis_vectors; ++j) {
            V1.item(0, j) = c[j];
        }
        proc_sampled_f_row[f_bv_max_global.proc].insert(f_bv_max_global.row);
        proc_f_row_to_tmp_fs_row[f_bv_max_global.proc][f_bv_max_global.row] = 0;
        num_samples_obtained++;
    }

    if (num_samples_obtained < num_samples)
    {
        Vector* A = new Vector(num_rows, f_basis->distributed());
        Vector* noM = new Vector(num_rows, f_basis->distributed());

        Matrix A0(num_basis_vectors - 1, num_basis_vectors - 1, false);
        Matrix V1_last_col(num_basis_vectors - 1, 1, false);
        Matrix tt(num_rows, num_basis_vectors - 1, f_basis->distributed());
        Matrix tt1(num_rows, num_basis_vectors - 1, f_basis->distributed());
        Matrix g1(tt.numRows(), tt.numColumns(), f_basis->distributed());
        Matrix GG(tt1.numRows(), tt1.numColumns(), f_basis->distributed());
        Vector ls_res_first_row(num_basis_vectors - 1, false);
        Vector nV(num_basis_vectors, false);

        int start_idx = 2;
        if (total_num_init_samples > 1)
        {
            start_idx += (total_num_init_samples - 1);
        }
        for (int i = start_idx; i <= num_samples; i++)
        {
            if (i <= num_basis_vectors)
            {
                A0.setSize(i - 1, i - 1);
                V1_last_col.setSize(i - 1, 1);
                for (int j = 0; j < num_samples_obtained; j++)
                {
                    for (int k = 0; k < i - 1; k++)
                    {
                        A0.item(j, k) = V1.item(j, k);
                    }
                    V1_last_col.item(j, 0) = V1.item(j, i - 1);
                }

                Matrix* atA0 = V1_last_col.transposeMult(A0);
                tt.setSize(num_rows, i - 1);
                tt1.setSize(num_rows, i - 1);

                double ata = 0.0;
                for (int j = 0; j < V1_last_col.numRows(); j++)
                {
                    for (int k = 0; k < V1_last_col.numColumns(); k++)
                    {
                        ata += (V1_last_col.item(j, k) * V1_last_col.item(j, k));
                    }
                }

                Matrix* lhs = A0.transposeMult(A0);
                lhs->inverse();
                Matrix* rhs = NULL;
                if (myid == 0)
                {
                    rhs = new Matrix(num_rows + atA0->numRows(), i - 1, f_basis->distributed());
                    for (int k = 0; k < rhs->numColumns(); k++)
                    {
                        rhs->item(0, k) = atA0->item(0, k);
                    }
                    for (int j = 1; j < rhs->numRows(); j++)
                    {
                        for (int k = 0; k < rhs->numColumns(); k++)
                        {
                            rhs->item(j, k) = Vo->item(j - 1, k);
                        }
                    }
                }
                else
                {
                    rhs = new Matrix(num_rows, i - 1, f_basis->distributed());
                    for (int j = 0; j < rhs->numRows(); j++)
                    {
                        for (int k = 0; k < rhs->numColumns(); k++)
                        {
                            rhs->item(j, k) = Vo->item(j, k);
                        }
                    }
                }

                Matrix* ls_res = rhs->mult(lhs);
                delete lhs;
                delete rhs;

                Matrix* c_T = NULL;
                if (myid == 0)
                {
                    c_T = new Matrix(ls_res->getData() + ls_res->numColumns(),
                                     ls_res->numRows() - 1, ls_res->numColumns(), f_basis->distributed(), true);
                }
                else
                {
                    c_T = new Matrix(ls_res->getData(),
                                     ls_res->numRows(), ls_res->numColumns(), f_basis->distributed(), true);
                }
                Matrix* Vo_first_i_columns = Vo->getFirstNColumns(i - 1);

                Vector* b = new Vector(num_rows, f_basis->distributed());
                for (int j = 0; j < Vo_first_i_columns->numRows(); j++)
                {
                    double tmp = 1.0;
                    for (int k = 0; k < Vo_first_i_columns->numColumns(); k++)
                    {
                        tmp += (Vo_first_i_columns->item(j, k) * c_T->item(j, k));
                    }
                    b->item(j) = tmp;
                }

                delete Vo_first_i_columns;

                for (int j = 0; j < num_rows; j++)
                {
                    for (int zz = 0; zz < i - 1; zz++)
                    {
                        tt.item(j, zz) = Vo->item(j, zz) * Vo->item(j, i - 1);
                    }
                }

                g1.setSize(tt.numRows(), tt.numColumns());
                for (int j = 0; j < g1.numRows(); j++)
                {
                    for (int k = 0; k < g1.numColumns(); k++)
                    {
                        g1.item(j, k) = atA0->item(0, k) + tt.item(j, k);
                    }
                }

                delete atA0;

                Vector* g3 = new Vector(num_rows, f_basis->distributed());
                for (int j = 0; j < c_T->numRows(); j++)
                {
                    double tmp = 0.0;
                    for (int k = 0; k < c_T->numColumns(); k++)
                    {
                        tmp += c_T->item(j, k) * g1.item(j, k);
                    }
                    g3->item(j) = tmp / b->item(j);
                }

                for (int j = 0; j < num_rows; j++)
                {
                    for (int zz = 0; zz < i - 1; zz++)
                    {
                        tt1.item(j, zz) = c_T->item(j, zz) * (Vo->item(j, i - 1) - g3->item(j));
                    }
                }

                delete c_T;
                delete g3;

                ls_res_first_row.setSize(ls_res->numColumns());
                if (myid == 0) {
                    for (int j = 0; j < ls_res->numColumns(); ++j) {
                        c[j] = ls_res->item(0, j);
                    }
                }
                MPI_Bcast(c.data(), ls_res->numColumns(), MPI_DOUBLE,
                          0, MPI_COMM_WORLD);
                for (int j = 0; j < ls_res->numColumns(); ++j) {
                    ls_res_first_row.item(j) = c[j];
                }
                GG.setSize(tt1.numRows(), tt1.numColumns());
                for (int j = 0; j < GG.numRows(); j++)
                {
                    for (int k = 0; k < GG.numColumns(); k++)
                    {
                        GG.item(j, k) = ls_res_first_row.item(k) + tt1.item(j, k);
                    }
                }

                delete ls_res;

                for (int j = 0; j < A->dim(); j++)
                {
                    double tmp = 0.0;
                    for (int k = 0; k < g1.numColumns(); k++)
                    {
                        tmp += g1.item(j, k) * GG.item(j, k);
                    }
                    A->item(j) = std::max(0.0, ata + (Vo->item(j, i - 1) * Vo->item(j,
                                                      i - 1)) - tmp);
                }

                nV.setSize(i);
                for (int j = 0; j < i; j++)
                {
                    nV.item(j) = 0.0;
                    for (int k = 0; k < num_samples_obtained; k++)
                    {
                        nV.item(j) += (V1.item(k, j) * V1.item(k, j));
                    }
                }

                for (int j = 0; j < noM->dim(); j++)
                {
                    noM->item(j) = 0.0;
                    for (int k = 0; k < i; k++)
                    {
                        noM->item(j) += std::log(nV.item(k) + (Vo->item(j, k) * Vo->item(j, k)));
                    }
                }

                for (int j = 0; j < A->dim(); j++)
                {
                    A->item(j) = std::log(fabs(A->item(j))) + std::log(b->item(j)) - noM->item(j);
                }

                delete b;
            }
            else
            {
                Matrix* curr_V1 = new Matrix(V1.getData(), num_samples_obtained,
                                             num_basis_vectors, false, true);
                Matrix* lhs = curr_V1->transposeMult(curr_V1);
                lhs->inverse();

                delete curr_V1;

                Matrix* ls_res = Vo->mult(lhs);
                delete lhs;

                nV.setSize(num_basis_vectors);
                for (int j = 0; j < num_basis_vectors; j++)
                {
                    nV.item(j) = 0.0;
                    for (int k = 0; k < num_samples_obtained; k++)
                    {
                        nV.item(j) += (V1.item(k, j) * V1.item(k, j));
                    }
                }

                for (int j = 0; j < noM->dim(); j++)
                {
                    noM->item(j) = 0.0;
                    for (int k = 0; k < num_basis_vectors; k++)
                    {
                        noM->item(j) += std::log(nV.item(k) + (Vo->item(j, k) * Vo->item(j, k)));
                    }
                }

                for (int j = 0; j < A->dim(); j++)
                {
                    double tmp = 0.0;
                    for (int k = 0; k < Vo->numColumns(); k++)
                    {
                        tmp += Vo->item(j, k) * ls_res->item(j, k);
                    }
                    A->item(j) = std::log(1 + tmp) - noM->item(j);
                }

                delete ls_res;
            }

            f_bv_max_local.row_val = -DBL_MAX;
            f_bv_max_local.proc = myid;
            for (int j = 0; j < num_rows; ++j) {
                if (proc_f_row_to_tmp_fs_row[f_bv_max_local.proc].find(j) ==
                        proc_f_row_to_tmp_fs_row[f_bv_max_local.proc].end())
                {
                    double f_bv_val = A->item(j);
                    if (f_bv_val > f_bv_max_local.row_val) {
                        f_bv_max_local.row_val = f_bv_val;
                        f_bv_max_local.row = j;
                    }
                }
            }
            MPI_Allreduce(&f_bv_max_local, &f_bv_max_global, 1,
                          MaxRowType, RowInfoOp, MPI_COMM_WORLD);
            // Now get the first sampled row of the basis
            if (f_bv_max_global.proc == myid) {
                for (int j = 0; j < num_basis_vectors; ++j) {
                    c[j] = Vo->item(f_bv_max_global.row, j);
                }
            }
            MPI_Bcast(c.data(), num_basis_vectors, MPI_DOUBLE,
                      f_bv_max_global.proc, MPI_COMM_WORLD);
            // Now add the first sampled row of the basis to tmp_fs.
            for (int j = 0; j < num_basis_vectors; ++j) {
                V1.item(num_samples_obtained, j) = c[j];
            }
            proc_sampled_f_row[f_bv_max_global.proc].insert(f_bv_max_global.row);
            proc_f_row_to_tmp_fs_row[f_bv_max_global.proc][f_bv_max_global.row] =
                num_samples_obtained;
            num_samples_obtained++;
        }

        delete A;
        delete noM;
    }

    // Fill f_sampled_row, and f_sampled_rows_per_proc.  Unscramble tmp_fs into
    // f_basis_sampled_inv.
    int idx = 0;
    for (int i = 0; i < num_procs; ++i) {
        set<int>& this_proc_sampled_f_row = proc_sampled_f_row[i];
        map<int, int>& this_proc_f_row_to_tmp_fs_row =
            proc_f_row_to_tmp_fs_row[i];
        f_sampled_rows_per_proc[i] = this_proc_sampled_f_row.size();
        for (set<int>::iterator j = this_proc_sampled_f_row.begin();
                j != this_proc_sampled_f_row.end(); ++j) {
            int this_f_row = *j;
            f_sampled_row[idx] = this_f_row;
            int tmp_fs_row = this_proc_f_row_to_tmp_fs_row[this_f_row];
            for (int col = 0; col < num_f_basis_cols; ++col) {
                f_basis_sampled_inv.item(idx, col) = V1.item(tmp_fs_row, col);
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

    if (qr_factorize)
    {
        delete Vo;
    }
    if (num_basis_vectors < f_basis->numColumns())
    {
        delete f_basis_truncated;
    }
}

}
