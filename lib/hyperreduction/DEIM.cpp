/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Interface to the DEIM algorithm to determine the rows of the
// rhs to be sampled for the interpolation of the rhs.

#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include "Utilities.h"
#include "mpi.h"
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <cfloat>
#include "DEIM.h"

using namespace std;

namespace CAROM {

void
DEIM(const Matrix* f_basis,
     int num_f_basis_vectors_used,
     std::vector<int>& f_sampled_row,
     std::vector<int>& f_sampled_rows_per_proc,
     Matrix& f_basis_sampled_inv,
     int myid,
     int num_procs,
     bool precond,
     Vector* K)
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
    int num_basis_vectors =
        std::min(num_f_basis_vectors_used, f_basis->numColumns());
    CAROM_VERIFY(num_basis_vectors == f_sampled_row.size());
    CAROM_VERIFY(num_basis_vectors == f_basis_sampled_inv.numRows()
                 && num_basis_vectors == f_basis_sampled_inv.numColumns());
    CAROM_VERIFY(!f_basis_sampled_inv.distributed());
    int basis_size = f_basis->numRows();
//(1-1)
    CAROM::Matrix* Kf = NULL;
    const Matrix* Vo = NULL;
    if(precond)
    {
        Kf = f_basis->row_normalize();
        Vo = Kf -> getFirstNColumns(num_basis_vectors);
    }else{
        Vo = f_basis;
    }
    int num_samples= num_basis_vectors;
    Vector Kii(num_samples,false);

    // The small matrix inverted by the algorithm.  We'll allocate the largest
    // matrix we'll need and set its size at each step in the algorithm.
    Matrix M(num_basis_vectors, num_basis_vectors, false);

    // Scratch space used throughout the algorithm.
//(1-2)
    double* c = new double [num_basis_vectors+1];

    vector<set<int> > proc_sampled_f_row(num_procs);
    vector<map<int, int> > proc_f_row_to_tmp_fs_row(num_procs);
    int num_f_basis_cols = f_basis_sampled_inv.numColumns();
    Matrix tmp_fs(f_basis_sampled_inv.numRows(),
                  num_f_basis_cols,
                  f_basis_sampled_inv.distributed());

    // Figure out the 1st sampled row of the RHS.
    RowInfo f_bv_max_local, f_bv_max_global;
    f_bv_max_local.row_val = -1.0;
    f_bv_max_local.proc = myid;
    for (int i = 0; i < basis_size; ++i) {
        double f_bv_val = fabs(f_basis->item(i, 0));
//(2-1)
	if(precond) f_bv_val = -Kf->item(i,num_basis_vectors);
        if(f_bv_val==-1.0) f_bv_val = -DBL_MAX;
        if (f_bv_val > f_bv_max_local.row_val) {
            f_bv_max_local.row_val = f_bv_val;
            f_bv_max_local.row = i;
        }
    }
    MPI_Allreduce(&f_bv_max_local, &f_bv_max_global, 1,
                  MaxRowType, RowInfoOp, MPI_COMM_WORLD);

    // Now get the first sampled row of the basis of the RHS.
    if (f_bv_max_global.proc == myid) {
	//(3-1)Vo
        for (int j = 0; j < num_basis_vectors; ++j) {
            c[j] = Vo->item(f_bv_max_global.row, j);
        }
        if(precond) c[num_basis_vectors] = Kf->item(f_bv_max_global.row, num_basis_vectors);
    }
    MPI_Bcast(c, num_basis_vectors+1, MPI_DOUBLE,
              f_bv_max_global.proc, MPI_COMM_WORLD);
    // Now add the first sampled row of the basis of the RHS to tmp_fs.
    for (int j = 0; j < num_basis_vectors; ++j) {
        tmp_fs.item(0, j) = c[j];
    }
//(3-2)
    if(precond) Kii.item(0) = c[num_basis_vectors];    
    proc_sampled_f_row[f_bv_max_global.proc].insert(f_bv_max_global.row);
    proc_f_row_to_tmp_fs_row[f_bv_max_global.proc][f_bv_max_global.row] = 0;

    // Now repeat the process for the other sampled rows of the basis of the
    // RHS.
    for (int i = 1; i < num_basis_vectors; ++i) {
        // If we currently know about S sampled rows of the basis of the RHS then
        // M contains the first S columns of those S sampled rows.
        M.setSize(i, i);
        for (int row = 0; row < i; ++row) {
            for (int col = 0; col < i; ++col) {
                M.item(row, col) = tmp_fs.item(row, col);
            }
        }

        // Invert M.
        M.inverse();

        // Now compute c, the inverse of M times the next column of the sampled
        // rows of the basis of the RHS.
        for (int minv_row = 0; minv_row < i; ++minv_row) {
            double tmp = 0.0;
            for (int minv_col = 0; minv_col < i; ++minv_col) {
                tmp += M.item(minv_row, minv_col)*tmp_fs.item(minv_col, i);
            }
            c[minv_row] = tmp;
        }

        // Now figure out the next sampled row of the basis of f.
        // Compute the first S basis vectors of the RHS times c and find the
        // row of this product have the greatest absolute value.  This is the
        // next sampled row of the basis of f.
        f_bv_max_local.row_val = -1.0;
        f_bv_max_local.proc = myid;
        for (int F_row = 0; F_row < basis_size; ++F_row) {
            double tmp = 0.0;
//(4-1)	    
            for (int F_col = 0; F_col < i; ++F_col) {
                tmp += Vo->item(F_row, F_col)*c[F_col];
            }
            double r_val = fabs(Vo->item(F_row, i) - tmp);
            if (r_val > f_bv_max_local.row_val) {
                f_bv_max_local.row_val = r_val;
                f_bv_max_local.row = F_row;
            }
        }
        MPI_Allreduce(&f_bv_max_local, &f_bv_max_global, 1,
                      MaxRowType, RowInfoOp, MPI_COMM_WORLD);

        // Now get the next sampled row of the basis of f.
        if (f_bv_max_global.proc == myid) {
//(4-2)	    
            for (int j = 0; j < num_basis_vectors; ++j) {
                c[j] = Vo->item(f_bv_max_global.row, j);
            }
	    if(precond) c[num_basis_vectors] = Kf->item(f_bv_max_global.row, num_basis_vectors);
        }
        MPI_Bcast(c, num_basis_vectors+1, MPI_DOUBLE,
                  f_bv_max_global.proc, MPI_COMM_WORLD);
        // Now add the ith sampled row of the basis of the RHS to tmp_fs.
        for (int j = 0; j < num_basis_vectors; ++j) {
            tmp_fs.item(i, j) = c[j];
        }
//(4-3)
	if(precond) Kii.item(i) = c[num_basis_vectors];
        proc_sampled_f_row[f_bv_max_global.proc].insert(f_bv_max_global.row);
        proc_f_row_to_tmp_fs_row[f_bv_max_global.proc][f_bv_max_global.row] = i;
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
                f_basis_sampled_inv.item(idx, col) = tmp_fs.item(tmp_fs_row, col);
            }
//(5-1)    
            if(precond) K->item(idx) = Kii.item(tmp_fs_row);	    
            ++idx;
        }
    }
//Temporary codes for printing results: should be erased 
    if(num_basis_vectors <= num_samples){
        for(int i=0; i<num_samples;i++){
                for(int j=0; j<num_basis_vectors; j++){
                //printf("%f\t",f_basis_sampled_inv.item(i,j));
                }
                //printf("\n");
        }
        Matrix* U = NULL;
        Matrix V(num_basis_vectors,num_basis_vectors,false);
        Vector sigma(num_basis_vectors, false);
        SerialSVD(&f_basis_sampled_inv, U, &sigma, &V);
        delete U;
        double sigma_end = 0.0;
        for(int j = num_basis_vectors-1; j>=0; j--){
                if(sigma.item(j) > 1e-8) {
                        sigma_end = sigma.item(j);
                        break;
                }
        }
        printf("conditionNum:%f -%f:%f\n", sigma.item(0)/sigma_end,sigma.item(0), sigma_end);
        printf("%d\t",f_basis_sampled_inv.numColumns());
        printf("%d\t",f_basis_sampled_inv.numRows());
    }
    if(precond){
      for(int i=0; i<num_samples; i++){
        printf("k%f\t",K->item(i));
      }
      printf("\n");
    }
    CAROM_ASSERT(num_samples == idx);

    // Now invert f_basis_sampled_inv.
    f_basis_sampled_inv.inverse();

    // Free the MPI_Datatype and MPI_Op.
    MPI_Type_free(&MaxRowType);
    MPI_Op_free(&RowInfoOp);
    if (precond)
    {
        delete Vo;
        delete Kf;
    }
    delete [] c;
}

}
