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
#define zgetrf CAROM_FC_GLOBAL(zgetrf, ZGETRF)
#define zgetri CAROM_FC_GLOBAL(zgetri, ZGETRI)

extern "C" {
    // Compute eigenvalue and eigenvectors of real non-symmetric matrix.
    void dgeev(char*, char*, int*, double*, int*, double*, double*, double*, int*, double*, int*, double*, int*, int*);

    // LU decomposition of a general matrix.
    void zgetrf(int*, int*, double*, int*, int*, int*);

    // Generate inverse of a matrix given its LU decomposition.
    void zgetri(int*, double*, int*, int*, double*, int*, int*);
}

using namespace std;

namespace CAROM {

DMD::DMD(int dim)
{
    // Get the rank of this process, and the number of processors.
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
    d_dim = dim;
}

bool DMD::takeSample(double* u_in)
{
    CAROM_VERIFY(u_in != 0);

    Vector sample(u_in, d_dim, true);
    d_snapshots.push_back(sample);
    return true;
}

void DMD::train(double energy_fraction)
{
    const Matrix* f_snapshots = getSnapshotMatrix();
    CAROM_VERIFY(f_snapshots->numColumns() > 1);
    CAROM_VERIFY(energy_fraction > 0 && energy_fraction <= 1);
    d_energy_fraction = energy_fraction;
    d_k = f_snapshots->numColumns() - 1;
    constructDMD(f_snapshots, d_rank, d_num_procs);

    delete f_snapshots;
}

void DMD::train(int k)
{
    const Matrix* f_snapshots = getSnapshotMatrix();
    CAROM_VERIFY(f_snapshots->numColumns() > 1);
    CAROM_VERIFY(k > 0 && k <= f_snapshots->numColumns() - 1);
    d_energy_fraction = -1;
    d_k = k;
    constructDMD(f_snapshots, d_rank, d_num_procs);

    delete f_snapshots;
}

void
DMD::constructDMD(const Matrix* f_snapshots,
                  int d_rank,
                  int d_num_procs)
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

    int *row_offset = new int[d_num_procs + 1];
    row_offset[d_num_procs] = f_snapshots_minus->numDistributedRows();
    row_offset[d_rank] = f_snapshots_minus->numRows();

    CAROM_VERIFY(MPI_Allgather(MPI_IN_PLACE,
                               1,
                               MPI_INT,
                               row_offset,
                               1,
                               MPI_INT,
                               MPI_COMM_WORLD) == MPI_SUCCESS);
    for (int i = d_num_procs - 1; i >= 0; i--) {
        row_offset[i] = row_offset[i + 1] - row_offset[i];
    }

    CAROM_VERIFY(row_offset[0] == 0);

    int d_blocksize = row_offset[d_num_procs] / d_num_procs;
    if (row_offset[d_num_procs] % d_num_procs != 0) d_blocksize += 1;

    SLPK_Matrix svd_input;

    // Calculate svd of snapshots_minus
    initialize_matrix(&svd_input, f_snapshots_minus->numColumns(),
                      f_snapshots_minus->numDistributedRows(),
                      1, d_num_procs, d_blocksize, d_blocksize);

    for (int rank = 0; rank < d_num_procs; ++rank)
    {
        scatter_block(&svd_input, 1, row_offset[rank] + 1,
                      f_snapshots_minus->getData(),
                      f_snapshots_minus->numColumns(),
                      row_offset[rank + 1] - row_offset[rank], rank);
    }

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
        int num_singular_vectors = std::min(f_snapshots_minus->numColumns(), f_snapshots_minus->numDistributedRows());
        for (int i = 0; i < num_singular_vectors; i++)
        {
            total_energy += d_factorizer->S[i];
        }
        double current_energy = 0.0;
        for (int i = 0; i < num_singular_vectors; i++)
        {
            current_energy += d_factorizer->S[i];
            if (current_energy / total_energy >= d_energy_fraction)
            {
                d_k = i + 1;
                break;
            }
        }
    }

    std::cout << "Using " << d_k << " basis vectors out of " << num_singular_vectors << "." << std::endl;

    // Allocate the appropriate matrices and gather their elements.
    Matrix* d_basis = new Matrix(f_snapshots->numRows(), d_k, f_snapshots->distributed());
    Matrix* d_S_inv = new Matrix(d_k, d_k, false);

    Matrix* d_basis_right = new Matrix(f_snapshots_minus->numColumns(), d_k, false);
    for (int d_rank = 0; d_rank < d_num_procs; ++d_rank) {
        // V is computed in the transposed order so no reordering necessary.
        gather_block(&d_basis->item(0, 0), d_factorizer->V,
                     1, row_offset[static_cast<unsigned>(d_rank)]+1,
                     d_k, row_offset[static_cast<unsigned>(d_rank) + 1] -
                     row_offset[static_cast<unsigned>(d_rank)],
                     d_rank);

        // gather_transposed_block does the same as gather_block, but transposes
        // it; here, it is used to go from column-major to row-major order.
        gather_transposed_block(&d_basis_right->item(0, 0), d_factorizer->U, 1, 1,
                                f_snapshots_minus->numColumns(), d_k, d_rank);
    }

    // Get inverse of singular values by multiplying by reciprocal.
    for (int i = 0; i < d_k; ++i)
    {
        d_S_inv->item(i, i) = 1 / d_factorizer->S[static_cast<unsigned>(i)];
    }

    // Calculate A_tilde = U_transpose * f_snapshots_plus * V * inv(S)
    Matrix* d_basis_mult_f_snapshots_plus = d_basis->transposeMult(f_snapshots_plus);
    Matrix* d_basis_mult_f_snapshots_plus_mult_d_basis_right = d_basis_mult_f_snapshots_plus->mult(d_basis_right);
    Matrix* A_tilde = d_basis_mult_f_snapshots_plus_mult_d_basis_right->mult(d_S_inv);

    // Calculate eigenvalues/eigenvectors of A_tilde
    // Call lapack routines to do the eigensolve.
    char jobvl, jobrl;

    // Calculate right eigenvectors only.
    jobvl = 'N';
    jobrl = 'V';

    int info;
    int lwork = std::max(d_k*d_k, 10*d_k);
    double* work = new double [lwork];
    double* e_real = new double [d_k];
    double* e_imaginary = new double [d_k];
    double* ev_l = NULL;
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
    dgeev(&jobvl, &jobrl, &d_k, A_tilde->getData(), &d_k, e_real, e_imaginary, ev_l, &d_k, ev_r->getData(), &d_k, work, &lwork, &info);

    // Eigenvalues now in a column major representation.  Put it
    // into row major order.
    for (int row = 0; row < d_k; ++row) {
        for (int col = row+1; col < d_k; ++col) {
            double tmp = ev_r->item(row, col);
            ev_r->item(row, col) = ev_r->item(col, row);
            ev_r->item(col, row) = tmp;
        }
    }

    Matrix* ev_real = new Matrix(d_k, d_k, false);
    Matrix* ev_imaginary = new Matrix(d_k, d_k, false);

    // Separate lapack eigenvector into real and imaginary parts
    for (int k = 0; k < d_k; ++k)
    {
        for (int row = 0; row < d_k; ++row) {
            ev_real->item(row, k) = ev_r->item(row, k);
        }
        if (e_imaginary[k] != 0)
        {
            for (int row = 0; row < d_k; ++row) {
                ev_real->item(row, k + 1) = ev_r->item(row, k);
                ev_imaginary->item(row, k) = ev_r->item(row, k + 1);
                ev_imaginary->item(row, k + 1) = -ev_r->item(row, k + 1);
            }

            // Skip the next eigenvalue since it'll be part of the complex
            // conjugate pair.
            ++k;
        }
    }

    for (int i = 0; i < d_k; i++)
    {
        d_eigs.push_back(std::complex<double>(e_real[i], e_imaginary[i]));
    }

    // Calculate phi
    Matrix* f_snapshots_plus_mult_d_basis_right = f_snapshots_plus->mult(d_basis_right);
    Matrix* f_snapshots_plus_mult_d_basis_right_mult_d_S_inv = f_snapshots_plus_mult_d_basis_right->mult(d_S_inv);
    d_phi_real = f_snapshots_plus_mult_d_basis_right_mult_d_S_inv->mult(ev_real);
    d_phi_imaginary = f_snapshots_plus_mult_d_basis_right_mult_d_S_inv->mult(ev_imaginary);

    Vector* init = new Vector(f_snapshots_minus->numRows(), true);
    for (int i = 0; i < init->dim(); i++)
    {
        init->item(i) = f_snapshots_minus->item(i, 0);
    }

    // Calculate pinv(d_phi) * initial_condition.
    projectInitialCondition(init);

    delete d_basis;
    delete d_basis_right;
    delete d_S_inv;
    delete d_basis_mult_f_snapshots_plus;
    delete d_basis_mult_f_snapshots_plus_mult_d_basis_right;
    delete A_tilde;
    delete f_snapshots_plus_mult_d_basis_right;
    delete f_snapshots_plus_mult_d_basis_right_mult_d_S_inv;
    delete f_snapshots_minus;
    delete f_snapshots_plus;
    delete ev_l;
    delete ev_r;
    delete ev_real;
    delete ev_imaginary;
    delete [] e_real;
    delete [] e_imaginary;
    delete init;
}


void
DMD::projectInitialCondition(const Vector* init)
{
    Matrix* d_phi_real_squared = d_phi_real->transposeMult(d_phi_real);
    Matrix* d_phi_real_squared_2 = d_phi_imaginary->transposeMult(d_phi_imaginary);
    *d_phi_real_squared -= *d_phi_real_squared_2;

    Matrix* d_phi_imaginary_squared = d_phi_real->transposeMult(d_phi_imaginary);
    Matrix* d_phi_imaginary_squared_2 = d_phi_imaginary->transposeMult(d_phi_real);
    *d_phi_imaginary_squared += *d_phi_imaginary_squared_2;

    double* inverse_input = new double[d_phi_real_squared->numRows() * d_phi_real_squared->numColumns() * 2];
    for (int i = 0; i < d_phi_real_squared->numRows(); i++)
    {
        int k = 0;
        for (int j = 0; j < d_phi_real_squared->numColumns() * 2; j++)
        {
            if (j % 2 == 0)
            {
                inverse_input[d_phi_real_squared->numColumns() * 2 * i + j] = d_phi_real_squared->item(i, k);
            }
            else
            {
                inverse_input[d_phi_imaginary_squared->numColumns() * 2 * i + j] = d_phi_imaginary_squared->item(i, k);
                k++;
            }
        }
    }

    // Call lapack routines to do the inversion.
    // Set up some stuff the lapack routines need.
    int info;
    int mtx_size = d_phi_real_squared->numColumns();
    int lwork = mtx_size*mtx_size;
    int* ipiv = new int [mtx_size];
    double* work = new double [lwork];

    // Now call lapack to do the inversion.
    zgetrf(&mtx_size, &mtx_size, inverse_input, &mtx_size, ipiv, &info);
    zgetri(&mtx_size, inverse_input, &mtx_size, ipiv, work, &lwork, &info);

    for (int i = 0; i < d_phi_real_squared->numRows(); i++)
    {
        int k = 0;
        for (int j = 0; j < d_phi_real_squared->numColumns() * 2; j++)
        {
            if (j % 2 == 0)
            {
                d_phi_real_squared->item(i, k) = inverse_input[d_phi_real_squared->numColumns() * 2 * i + j];
            }
            else
            {
                d_phi_imaginary_squared->item(i, k) = inverse_input[d_phi_imaginary_squared->numColumns() * 2 * i + j];
                k++;
            }
        }
    }

    Vector* rhs_real = d_phi_real->transposeMult(init);
    Vector* rhs_imaginary = d_phi_imaginary->transposeMult(init);

    Vector* d_projected_init_real_1 = d_phi_real_squared->mult( rhs_real);
    Vector* d_projected_init_real_2 = d_phi_imaginary_squared->mult(rhs_imaginary);
    d_projected_init_real = d_projected_init_real_1->minus(d_projected_init_real_2);

    Vector* d_projected_init_imaginary_1 = d_phi_real_squared->mult(rhs_imaginary);
    Vector* d_projected_init_imaginary_2 = d_phi_imaginary_squared->mult(rhs_real);
    d_projected_init_imaginary =  d_projected_init_imaginary_1->plus(d_projected_init_imaginary_2);

    delete d_phi_real_squared;
    delete d_phi_real_squared_2;
    delete  d_projected_init_real_1;
    delete  d_projected_init_real_2;
    delete d_phi_imaginary_squared;
    delete d_phi_imaginary_squared_2;
    delete  d_projected_init_imaginary_1;
    delete  d_projected_init_imaginary_2;
    delete rhs_real;
    delete rhs_imaginary;

    delete [] inverse_input;
    delete [] ipiv;
    delete [] work;
}

Vector*
DMD::predict(double t,
             double dt)
{
    const std::pair<Vector*, Vector*> d_projected_init_pair(d_projected_init_real, d_projected_init_imaginary);
    return predict(d_projected_init_pair, t, dt);
}

Vector*
DMD::predict(const std::pair<Vector*, Vector*> init,
             double t,
             double dt)
{
    std::pair<Matrix*, Matrix*> d_phi_pair = phiMultEigs(t, dt);
    Matrix* d_phi_mult_eigs_real = d_phi_pair.first;
    Matrix* d_phi_mult_eigs_imaginary = d_phi_pair.second;

    Vector* d_predicted_state_real_1 = d_phi_mult_eigs_real->mult(init.first);
    Vector* d_predicted_state_real_2 = d_phi_mult_eigs_imaginary->mult(init.second);
    Vector* d_predicted_state_real = d_predicted_state_real_1->minus(d_predicted_state_real_2);

    delete d_phi_mult_eigs_real;
    delete d_phi_mult_eigs_imaginary;
    delete d_predicted_state_real_1;
    delete d_predicted_state_real_2;

    return d_predicted_state_real;
}

std::pair<Matrix*, Matrix*>
DMD::phiMultEigs(double t,
                 double dt)
{
    Matrix* d_eigs_exp_real = new Matrix(d_k, d_k, false);
    Matrix* d_eigs_exp_imaginary = new Matrix(d_k, d_k, false);

    double exp_power = t / dt;

    for (int i = 0; i < d_k; i++)
    {
        std::complex<double> eig_exp = std::pow(d_eigs[i], exp_power);
        d_eigs_exp_real->item(i, i) = std::real(eig_exp);
        d_eigs_exp_imaginary->item(i, i) = std::imag(eig_exp);
    }

    Matrix* d_phi_mult_eigs_real = d_phi_real->mult(d_eigs_exp_real);
    Matrix* d_phi_mult_eigs_real_2 = d_phi_imaginary->mult(d_eigs_exp_imaginary);
    *d_phi_mult_eigs_real -= *d_phi_mult_eigs_real_2;
    Matrix* d_phi_mult_eigs_imaginary = d_phi_real->mult(d_eigs_exp_imaginary);
    Matrix* d_phi_mult_eigs_imaginary_2 = d_phi_imaginary->mult(d_eigs_exp_real);
    *d_phi_mult_eigs_imaginary += *d_phi_mult_eigs_imaginary_2;

    delete d_eigs_exp_real;
    delete d_eigs_exp_imaginary;
    delete d_phi_mult_eigs_real_2;
    delete d_phi_mult_eigs_imaginary_2;

    return std::pair<Matrix*,Matrix*>(d_phi_mult_eigs_real, d_phi_mult_eigs_imaginary);
}

const Matrix*
DMD::getSnapshotMatrix()
{
    CAROM_VERIFY(d_snapshots.size() > 0);
    CAROM_VERIFY(d_snapshots[0].dim() > 0);
    for (int i = 0 ; i < d_snapshots.size() - 1; i++)
    {
        CAROM_VERIFY(d_snapshots[i].dim() == d_snapshots[i + 1].dim());
        CAROM_VERIFY(d_snapshots[i].distributed() == d_snapshots[i + 1].distributed());
    }

    Matrix* snapshot_mat = new Matrix(d_snapshots[0].dim(), d_snapshots.size(), d_snapshots[0].distributed());

    // This might be slow since we're accessing a different memory address each time.
    for (int i = 0; i < d_snapshots[0].dim(); i++)
    {
        for (int j = 0; j < d_snapshots.size(); j++)
        {
            snapshot_mat->item(i, j) = d_snapshots[j].item(i);
        }
    }

    return snapshot_mat;
}

}
