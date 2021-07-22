/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Interface to the DMD algorithm.

#include "DMD.h"

#include "Matrix.h"
#include "Vector.h"
#include "scalapack_wrapper.h"
#include "mpi.h"

/* Use C++11 built-in shared pointers if available; else fallback to Boost. */
#if __cplusplus >= 201103L
#include <memory>
#else
#include <boost/shared_ptr.hpp>
#endif

#include "scalapack_wrapper.h"

/* Use automatically detected Fortran name-mangling scheme */
#define dgeev CAROM_FC_GLOBAL(dgeev, DGEEV)

extern "C" {
    // Compiute eigenvalue and eigenvectors of real non-symmetric matrix.
    void dgeev(char*, char*, int*, double*, int*, double*, double*, double*, int*, double*, int*, double*, int*, int*);
}

using namespace std;

namespace CAROM {

    DMD::DMD(const Matrix* f_snapshots,
             double energy_fraction,
             int rank,
             int num_procs)
    {
        CAROM_VERIFY(f_snapshots->numColumns() > 1);
        CAROM_VERIFY(energy_fraction > 0 && energy_fraction <= 1);
        d_energy_fraction = energy_fraction;
        d_k = f_snapshots->numColumns() - 1;
        constructDMD(f_snapshots, rank, num_procs);
    }

    DMD::DMD(const Matrix* f_snapshots,
             int k,
             int rank,
             int num_procs)
    {
        CAROM_VERIFY(f_snapshots->numColumns() > 1);
        CAROM_VERIFY(k > 0 && k <= f_snapshots->numColumns() - 1);
        d_energy_fraction = -1;
        d_k = k;
        constructDMD(f_snapshots, rank, num_procs);
    }

    void
    DMD::constructDMD(const Matrix* f_snapshots,
                      int rank,
                      int num_procs)
    {
        Matrix* f_snapshots_minus = new Matrix(f_snapshots->numRows(),
            f_snapshots->numColumns() - 1, f_snapshots->distributed());
        Matrix* f_snapshots_plus = new Matrix(f_snapshots->numRows(),
            f_snapshots->numColumns() - 1, f_snapshots->distributed());


        // Break up snapshots into snapshots_minus and snapshots_plus
        // snapshots_minus = all columns of snapshots except last
        // snapshots_plus = all columns of snapshots except first
        for (int i = 0; i < f_snapshots->numRows(); i++)
        {
            f_snapshots_minus->item(i, 0) = f_snapshots->item(i, 0);
            for (int j = 1; j < f_snapshots->numColumns() - 1; j++)
            {
                f_snapshots_minus->item(i, j) = f_snapshots->item(i, j);
                f_snapshots_plus->item(i, j - 1) = f_snapshots->item(i, j);
            }
            f_snapshots_plus->item(i, f_snapshots->numColumns() - 2) =
                f_snapshots->item(i, f_snapshots->numColumns() - 1);
        }

        int *row_offset = new int[num_procs + 1];
        row_offset[num_procs] = f_snapshots_minus->numDistributedRows();
        row_offset[rank] = f_snapshots_minus->numRows();

        CAROM_VERIFY(MPI_Allgather(MPI_IN_PLACE,
                                   1,
                                   MPI_INT,
                                   row_offset,
                                   1,
                                   MPI_INT,
                                   MPI_COMM_WORLD) == MPI_SUCCESS);
        for (int i = num_procs - 1; i >= 0; i--) {
            row_offset[i] = row_offset[i + 1] - row_offset[i];
        }

        CAROM_VERIFY(row_offset[0] == 0);

        int d_blocksize = row_offset[num_procs] / num_procs;
        if (row_offset[num_procs] % num_procs != 0) d_blocksize += 1;

        SLPK_Matrix svd_input;

        // Calculate svd of snapshots_minus
        initialize_matrix(&svd_input, f_snapshots_minus->numColumns(),
                          f_snapshots_minus->numDistributedRows(),
                          1, num_procs, d_blocksize, d_blocksize);
        scatter_block(&svd_input, 1, 1,
                      f_snapshots_minus->getData(),
                      f_snapshots_minus->numColumns(),
                      f_snapshots_minus->numDistributedRows(),0);

        std::unique_ptr<SVDManager> d_factorizer(new SVDManager);

        // This block does the actual ScaLAPACK call to do the factorization.
        svd_init(d_factorizer.get(), &svd_input);

        d_factorizer->dov = 1;

        factorize(d_factorizer.get());
        free_matrix_data(&svd_input);

        // Compute how many basis vectors we will actually use.
        if (d_energy_fraction != -1)
        {
            double total_energy = 0.0;
            for (int i = 0; i < f_snapshots->numColumns(); i++)
            {
                total_energy += d_factorizer->S[i];
            }
            double current_energy = 0.0;
            for (int i = 0; i < f_snapshots->numColumns(); i++)
            {
                current_energy += d_factorizer->S[i];
                if (current_energy / total_energy >= d_energy_fraction)
                {
                    d_k = i;
                    break;
                }
            }
        }

        // Allocate the appropriate matrices and gather their elements.
        Matrix* d_basis = new Matrix(f_snapshots->numRows(), d_k, true);
        Vector* d_S_inv = new Vector(d_k, false);

        Matrix* d_basis_right = new Matrix(d_k, d_k, false);
        for (int rank = 0; rank < num_procs; ++rank) {
            // gather_transposed_block does the same as gather_block, but transposes
            // it; here, it is used to go from column-major to row-major order.
            gather_transposed_block(&d_basis->item(0, 0), d_factorizer->U,
                                    row_offset[static_cast<unsigned>(rank)]+1,
                                    1, row_offset[static_cast<unsigned>(rank) + 1] -
                                    row_offset[static_cast<unsigned>(rank)],
                                    d_k, rank);
            // V is computed in the transposed order so no reordering necessary.
            gather_block(&d_basis_right->item(0, 0), d_factorizer->V, 1, 1,
                         d_k, d_k, rank);
        }

        // Get inverse of singular values by multiplying by reciprocal.
        for (int i = 0; i < d_k; ++i)
        {
            d_S_inv->item(i) = 1 / d_factorizer->S[static_cast<unsigned>(i)];
        }

        // Calculate A_tilde = U_transpose * f_snapshots_plus * V * inv(S)
        Matrix* A_tilde = d_basis->transposeMult(f_snapshots_plus);
        A_tilde->mult(d_basis_right);
        A_tilde->mult(d_S_inv);

        // Calculate eigenvalues/eigenvectors of A_tilde
        // Call lapack routines to do the eigensolve.
        char jobvl, jobrl;

        // Calculate right eigenvectors only.
        jobvl = 'N';
        jobrl = 'V';

        int info;
        int lwork = d_k * d_k;
        double* work = new double [lwork];
        double* e_real = new double [d_k];
        double* e_imaginary = new double [d_k];
        Matrix* ev_l = new Matrix(d_k, d_k, false);
        Matrix* ev_r = new Matrix(d_k, d_k, false);

        // A_tilde now in a row major representation.  Put it
        // into column major order.
        for (int row = 0; row < d_k; ++row) {
            for (int col = row+1; col < d_k; ++col) {
                double tmp = A_tilde->item(row, col);
                A_tilde->item(row, col) = A_tilde->item(col, row);
                A_tilde->item(col, row) = tmp;
            }
        }

        // Now call lapack to do the eigensolve.
        dgeev(&jobvl, &jobrl, &d_k, A_tilde->getData(), &d_k, e_real, e_imaginary, ev_l->getData(), &d_k, ev_r->getData(), &d_k, work, &lwork, &info);

        // Eigenvalues now in a column major representation.  Put it
        // into row major order.
        for (int row = 0; row < d_k; ++row) {
            for (int col = row+1; col < d_k; ++col) {
                double tmp = ev_r->item(row, col);
                ev_r->item(row, col) = ev_r->item(col, row);
                ev_r->item(col, row) = tmp;
            }
        }

        delete [] e_real;
        delete [] e_imaginary;

        // Calculate phi
        d_phi = f_snapshots_plus->mult(d_basis_right);
        d_phi->mult(d_S_inv);
        d_phi->mult(ev_r);

        delete d_basis;
        delete d_basis_right;
        delete d_S_inv;
        delete A_tilde;
        delete f_snapshots_minus;
        delete f_snapshots_plus;
        delete ev_l;
        delete ev_r;
    }

}
