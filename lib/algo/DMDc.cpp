/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Implementation of the DMDc algorithm.

#include "DMD.h"
#include "DMDc.h"

#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include "linalg/scalapack_wrapper.h"
#include "utils/Utilities.h"
#include "utils/CSVDatabase.h"
#include "utils/HDFDatabase.h"
#include "mpi.h"

#include <cstring>

/* Use C++11 built-in shared pointers if available; else fallback to Boost. */
#if __cplusplus >= 201103L
#include <memory>
#else
#include <boost/shared_ptr.hpp>
#endif

/* Use automatically detected Fortran name-mangling scheme */
#define zgetrf CAROM_FC_GLOBAL(zgetrf, ZGETRF)
#define zgetri CAROM_FC_GLOBAL(zgetri, ZGETRI)

extern "C" {
    // LU decomposition of a general matrix.
    void zgetrf(int*, int*, double*, int*, int*, int*);

    // Generate inverse of a matrix given its LU decomposition.
    void zgetri(int*, double*, int*, int*, double*, int*, int*);
}

namespace CAROM {

DMDc::DMDc(int dim, int dim_c, std::shared_ptr<Vector> state_offset)
{
    CAROM_VERIFY(dim > 0);
    CAROM_VERIFY(dim_c > 0);

    // Get the rank of this process, and the number of processors.
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
    d_dim = dim;
    d_dim_c = dim_c;
    d_trained = false;
    d_init_projected = false;
    setOffset(state_offset);
}

DMDc::DMDc(int dim, int dim_c, double dt, std::shared_ptr<Vector> state_offset)
{
    CAROM_VERIFY(dim > 0);
    CAROM_VERIFY(dim_c > 0);
    CAROM_VERIFY(dt > 0.0);

    // Get the rank of this process, and the number of processors.
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
    d_dim = dim;
    d_dim_c = dim_c;
    d_dt = dt;
    d_trained = false;
    d_init_projected = false;
    setOffset(state_offset);
}

DMDc::DMDc(std::string base_file_name)
{
    // Get the rank of this process, and the number of processors.
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
    d_trained = true;
    d_init_projected = true;

    load(base_file_name);
}

DMDc::DMDc(std::vector<std::complex<double>> & eigs,
           std::shared_ptr<Matrix> & phi_real,
           std::shared_ptr<Matrix> & phi_imaginary,
           std::shared_ptr<Matrix> & B_tilde, int k, double dt, double t_offset,
           std::shared_ptr<Vector> & state_offset, std::shared_ptr<Matrix> & basis)
{
    // Get the rank of this process, and the number of processors.
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
    d_trained = true;
    d_init_projected = false;

    d_eigs = eigs;
    d_phi_real = phi_real;
    d_phi_imaginary = phi_imaginary;
    d_B_tilde = B_tilde;
    d_k = k;
    d_dt = dt;
    d_t_offset = t_offset;
    d_basis = basis;
    setOffset(state_offset);
}

void DMDc::setOffset(std::shared_ptr<Vector> & offset_vector)
{
    d_state_offset = offset_vector;
}

void DMDc::takeSample(double* u_in, double t, double* f_in, bool last_step)
{
    CAROM_VERIFY(u_in != 0);
    CAROM_VERIFY(t >= 0.0);
    Vector* sample = new Vector(u_in, d_dim, true);

    double orig_t = t;
    if (d_snapshots.empty())
    {
        d_t_offset = t;
        t = 0.0;
    }
    else
    {
        t -= d_t_offset;
    }

    // Erase any snapshots taken at the same or later time
    while (!d_sampled_times.empty() && d_sampled_times.back() >= t)
    {
        if (d_rank == 0) std::cout << "Removing existing snapshot at time: " <<
                                       d_t_offset + d_sampled_times.back() << std::endl;
        d_snapshots.pop_back();
        d_controls.pop_back();
        d_sampled_times.pop_back();
    }

    if (d_snapshots.empty())
    {
        d_t_offset = orig_t;
        t = 0.0;
    }
    else
    {
        CAROM_VERIFY(d_sampled_times.back() < t);
    }
    d_snapshots.push_back(std::shared_ptr<Vector>(sample));

    if (!last_step)
    {
        Vector* control = new Vector(f_in, d_dim_c, false);
        d_controls.push_back(std::shared_ptr<Vector>(control));
    }

    d_sampled_times.push_back(t);
}

void DMDc::train(double energy_fraction, const Matrix* B)
{
    std::unique_ptr<const Matrix> f_snapshots = getSnapshotMatrix();
    std::unique_ptr<const Matrix> f_controls = createSnapshotMatrix(d_controls);
    CAROM_VERIFY(f_snapshots->numColumns() > 1);
    CAROM_VERIFY(f_controls->numColumns() == f_snapshots->numColumns() - 1);
    CAROM_VERIFY(energy_fraction > 0 && energy_fraction <= 1);
    d_energy_fraction = energy_fraction;
    constructDMDc(*f_snapshots, *f_controls, d_rank, d_num_procs, B);
}

void DMDc::train(int k, const Matrix* B)
{
    std::unique_ptr<const Matrix> f_snapshots = getSnapshotMatrix();
    std::unique_ptr<const Matrix> f_controls = createSnapshotMatrix(d_controls);
    CAROM_VERIFY(f_snapshots->numColumns() > 1);
    CAROM_VERIFY(f_controls->numColumns() == f_snapshots->numColumns() - 1);
    CAROM_VERIFY(k > 0 && k <= f_snapshots->numColumns() - 1);
    d_energy_fraction = -1.0;
    d_k = k;
    constructDMDc(*f_snapshots, *f_controls, d_rank, d_num_procs, B);
}

std::pair<Matrix*, Matrix*>
DMDc::computeDMDcSnapshotPair(const Matrix & snapshots, const Matrix & controls,
                              const Matrix* B)
{
    CAROM_VERIFY(snapshots.numColumns() > 1);
    CAROM_VERIFY(controls.numColumns() == snapshots.numColumns() - 1);

    // TODO: Making two copies of the snapshot matrix has a lot of overhead.
    //       We need to figure out a way to do submatrix multiplication and to
    //       reimplement this algorithm using one snapshot matrix.
    int input_control_dim = (B == NULL) ? controls.numRows() : 0;
    Matrix* f_snapshots_in = new Matrix(snapshots.numRows() + input_control_dim,
                                        snapshots.numColumns() - 1, snapshots.distributed());
    Matrix* f_snapshots_out = new Matrix(snapshots.numRows(),
                                         snapshots.numColumns() - 1, snapshots.distributed());

    // Break up snapshots into snapshots_in and snapshots_out
    // snapshots_in = all columns of snapshots except last
    // snapshots_out = all columns of snapshots except first
    for (int i = 0; i < snapshots.numRows(); i++)
    {
        for (int j = 0; j < snapshots.numColumns() - 1; j++)
        {
            f_snapshots_in->item(i, j) = snapshots(i, j);
            f_snapshots_out->item(i, j) = snapshots(i, j + 1);
            if (d_state_offset)
            {
                f_snapshots_in->item(i, j) -= d_state_offset->item(i);
                f_snapshots_out->item(i, j) -= d_state_offset->item(i);
            }
        }
    }

    for (int i = 0; i < controls.numRows(); i++)
    {
        if (B == NULL)
        {
            for (int j = 0; j < snapshots.numColumns() - 1; j++)
            {
                f_snapshots_in->item(snapshots.numRows() + i, j) = controls.item(i, j);
            }
        }
        else
        {
            std::unique_ptr<Matrix> Bf = B->mult(controls);
            *f_snapshots_out -= *Bf;
        }
    }

    return std::pair<Matrix*,Matrix*>(f_snapshots_in, f_snapshots_out);
}

void
DMDc::constructDMDc(const Matrix & f_snapshots,
                    const Matrix & f_controls,
                    int d_rank,
                    int d_num_procs,
                    const Matrix* B)
{
    std::pair<Matrix*, Matrix*> f_snapshot_pair = computeDMDcSnapshotPair(
                f_snapshots, f_controls, B);
    Matrix* f_snapshots_in = f_snapshot_pair.first;
    Matrix* f_snapshots_out = f_snapshot_pair.second;

    int *row_offset = new int[d_num_procs + 1];
    row_offset[d_num_procs] = f_snapshots_in->numDistributedRows();
    row_offset[d_rank] = f_snapshots_in->numRows();

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

    // Calculate svd of snapshots_in
    initialize_matrix(&svd_input, f_snapshots_in->numColumns(),
                      f_snapshots_in->numDistributedRows(),
                      1, d_num_procs, d_blocksize, d_blocksize);

    for (int rank = 0; rank < d_num_procs; ++rank)
    {
        scatter_block(&svd_input, 1, row_offset[rank] + 1,
                      f_snapshots_in->getData(),
                      f_snapshots_in->numColumns(),
                      row_offset[rank + 1] - row_offset[rank], rank);
    }

    std::unique_ptr<SVDManager> d_factorizer_in(new SVDManager);

    // This block does the actual ScaLAPACK call to do the factorization.
    svd_init(d_factorizer_in.get(), &svd_input);
    d_factorizer_in->dov = 1;
    factorize(d_factorizer_in.get());
    free_matrix_data(&svd_input);

    // Compute how many basis vectors we will actually use.
    d_num_singular_vectors = std::min(f_snapshots_in->numColumns(),
                                      f_snapshots_in->numDistributedRows());
    for (int i = 0; i < d_num_singular_vectors; i++)
    {
        d_sv.push_back(d_factorizer_in->S[i]);
    }

    int d_k_in = d_k;
    if (d_energy_fraction != -1.0)
    {
        d_k_in = d_num_singular_vectors;
        if (d_energy_fraction < 1.0)
        {
            double total_energy = 0.0;
            for (int i = 0; i < d_num_singular_vectors; i++)
            {
                total_energy += d_factorizer_in->S[i];
            }
            double current_energy = 0.0;
            for (int i = 0; i < d_num_singular_vectors; i++)
            {
                current_energy += d_factorizer_in->S[i];
                if (current_energy / total_energy >= d_energy_fraction)
                {
                    d_k_in = i + 1;
                    break;
                }
            }
        }
    }

    if (d_rank == 0) std::cout << "Using " << d_k_in << " basis vectors out of " <<
                                   d_num_singular_vectors << " for input." << std::endl;

    // Allocate the appropriate matrices and gather their elements.
    std::unique_ptr<Matrix> d_basis_in(new Matrix(f_snapshots_in->numRows(), d_k_in,
                                       f_snapshots_in->distributed()));
    Matrix* d_S_inv = new Matrix(d_k_in, d_k_in, false);
    Matrix* d_basis_right = new Matrix(f_snapshots_in->numColumns(), d_k_in, false);

    for (int d_rank = 0; d_rank < d_num_procs; ++d_rank) {
        // V is computed in the transposed order so no reordering necessary.
        gather_block(&d_basis_in->item(0, 0), d_factorizer_in->V,
                     1, row_offset[static_cast<unsigned>(d_rank)]+1,
                     d_k_in, row_offset[static_cast<unsigned>(d_rank) + 1] -
                     row_offset[static_cast<unsigned>(d_rank)],
                     d_rank);

        // gather_transposed_block does the same as gather_block, but transposes
        // it; here, it is used to go from column-major to row-major order.
        gather_transposed_block(&d_basis_right->item(0, 0), d_factorizer_in->U, 1, 1,
                                f_snapshots_in->numColumns(), d_k_in, d_rank);
    }

    // Get inverse of singular values by multiplying by reciprocal.
    for (int i = 0; i < d_k_in; ++i)
    {
        d_S_inv->item(i, i) = 1 / d_factorizer_in->S[static_cast<unsigned>(i)];
    }

    if (B == NULL)
    {
        // SVD on outputs
        row_offset[d_num_procs] = f_snapshots_out->numDistributedRows();
        row_offset[d_rank] = f_snapshots_out->numRows();

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

        SLPK_Matrix svd_output;

        // Calculate svd of snapshots_out
        initialize_matrix(&svd_output, f_snapshots_out->numColumns(),
                          f_snapshots_out->numDistributedRows(),
                          1, d_num_procs, d_blocksize, d_blocksize);

        for (int rank = 0; rank < d_num_procs; ++rank)
        {
            scatter_block(&svd_output, 1, row_offset[rank] + 1,
                          f_snapshots_out->getData(),
                          f_snapshots_out->numColumns(),
                          row_offset[rank + 1] - row_offset[rank], rank);
        }

        std::unique_ptr<SVDManager> d_factorizer_out(new SVDManager);

        // This block does the actual ScaLAPACK call to do the factorization.
        svd_init(d_factorizer_out.get(), &svd_output);
        d_factorizer_out->dov = 1;
        factorize(d_factorizer_out.get());
        free_matrix_data(&svd_output);

        // Compute how many basis vectors we will actually use.
        d_num_singular_vectors = std::min(f_snapshots_out->numColumns(),
                                          f_snapshots_out->numDistributedRows());
        for (int i = 0; i < d_num_singular_vectors; i++)
        {
            d_sv.push_back(d_factorizer_out->S[i]);
        }

        if (d_energy_fraction != -1.0)
        {
            d_k = d_num_singular_vectors;
            if (d_energy_fraction < 1.0)
            {
                double total_energy = 0.0;
                for (int i = 0; i < d_num_singular_vectors; i++)
                {
                    total_energy += d_factorizer_out->S[i];
                }
                double current_energy = 0.0;
                for (int i = 0; i < d_num_singular_vectors; i++)
                {
                    current_energy += d_factorizer_out->S[i];
                    if (current_energy / total_energy >= d_energy_fraction)
                    {
                        d_k = i + 1;
                        break;
                    }
                }
            }
        }

        if (d_rank == 0) std::cout << "Using " << d_k << " basis vectors out of " <<
                                       d_num_singular_vectors << " for output." << std::endl;

        // Allocate the appropriate matrices and gather their elements.
        d_basis.reset(new Matrix(f_snapshots_out->numRows(), d_k,
                                 f_snapshots_out->distributed()));
        for (int d_rank = 0; d_rank < d_num_procs; ++d_rank) {
            // V is computed in the transposed order so no reordering necessary.
            gather_block(&d_basis->item(0, 0), d_factorizer_out->V,
                         1, row_offset[static_cast<unsigned>(d_rank)]+1,
                         d_k, row_offset[static_cast<unsigned>(d_rank) + 1] -
                         row_offset[static_cast<unsigned>(d_rank)],
                         d_rank);
        }

        free_matrix_data(d_factorizer_out->U);
        free_matrix_data(d_factorizer_out->V);
        free(d_factorizer_out->U);
        free(d_factorizer_out->V);
        free(d_factorizer_out->S);

        release_context(&svd_output);
    }
    else
    {
        d_basis.reset(d_basis_in.release());
        d_k = d_k_in;
    }

    delete[] row_offset;

    // Calculate A_tilde and B_tilde
    std::unique_ptr<Matrix> d_basis_mult_f_snapshots_out = d_basis->transposeMult(
                *f_snapshots_out);
    std::unique_ptr<Matrix> d_basis_mult_f_snapshots_out_mult_d_basis_right =
        d_basis_mult_f_snapshots_out->mult(*d_basis_right);
    std::shared_ptr<Matrix> d_A_tilde_orig =
        d_basis_mult_f_snapshots_out_mult_d_basis_right->mult(
            *d_S_inv);

    if (B == NULL)
    {
        Matrix* d_basis_in_state = new Matrix(f_snapshots.numRows(),
                                              d_k_in, f_snapshots.distributed());
        Matrix* d_basis_in_control_transpose = new Matrix(d_k_in, f_controls.numRows(),
                false);
        for (int j = 0; j < d_k_in; j++)
        {
            for (int i = 0; i < f_snapshots.numRows(); i++)
            {
                d_basis_in_state->item(i, j) = d_basis_in->item(i, j);
            }
            for (int i = 0; i < f_controls.numRows(); i++)
            {
                d_basis_in_control_transpose->item(j, i) =
                    d_basis_in->item(f_snapshots.numRows() + i, j);
            }
        }
        std::unique_ptr<Matrix> d_basis_state_rot = d_basis_in_state->transposeMult(
                    *d_basis);
        d_A_tilde = d_A_tilde_orig->mult(*d_basis_state_rot);
        d_B_tilde = d_A_tilde_orig->mult(*d_basis_in_control_transpose);
        delete d_basis_in_state;
        delete d_basis_in_control_transpose;
    }
    else
    {
        d_A_tilde = d_A_tilde_orig;
        d_B_tilde = d_basis->transposeMult(*B);
    }

    // Calculate the right eigenvalues/eigenvectors of A_tilde
    ComplexEigenPair eigenpair = NonSymmetricRightEigenSolve(*d_A_tilde);
    d_eigs = eigenpair.eigs;

    d_phi_real = d_basis->mult(*eigenpair.ev_real);
    d_phi_imaginary = d_basis->mult(*eigenpair.ev_imaginary);

    Vector init(d_basis->numRows(), true);
    for (int i = 0; i < init.dim(); i++)
    {
        init(i) = f_snapshots_in->item(i, 0);
    }

    // Calculate the projection initial_condition onto column space of d_basis.
    project(init, f_controls);

    d_trained = true;

    delete d_basis_right;
    delete d_S_inv;
    delete f_snapshots_in;
    delete f_snapshots_out;
    delete eigenpair.ev_real;
    delete eigenpair.ev_imaginary;

    free_matrix_data(d_factorizer_in->U);
    free_matrix_data(d_factorizer_in->V);
    free(d_factorizer_in->U);
    free(d_factorizer_in->V);
    free(d_factorizer_in->S);

    release_context(&svd_input);
}

void
DMDc::project(const Vector & init, const Matrix & controls, double t_offset)
{
    std::shared_ptr<Matrix> d_phi_real_squared = d_phi_real->transposeMult(
                *d_phi_real);
    std::unique_ptr<Matrix> d_phi_real_squared_2 = d_phi_imaginary->transposeMult(
                *d_phi_imaginary);
    *d_phi_real_squared += *d_phi_real_squared_2;

    std::shared_ptr<Matrix> d_phi_imaginary_squared = d_phi_real->transposeMult(
                *d_phi_imaginary);
    std::unique_ptr<Matrix> d_phi_imaginary_squared_2 =
        d_phi_imaginary->transposeMult(*d_phi_real);
    *d_phi_imaginary_squared -= *d_phi_imaginary_squared_2;

    const int dprs_row = d_phi_real_squared->numRows();
    const int dprs_col = d_phi_real_squared->numColumns();
    double* inverse_input = new double[dprs_row * dprs_col * 2];
    for (int i = 0; i < d_phi_real_squared->numRows(); i++)
    {
        int k = 0;
        for (int j = 0; j < d_phi_real_squared->numColumns() * 2; j++)
        {
            if (j % 2 == 0)
            {
                inverse_input[d_phi_real_squared->numColumns() * 2 * i + j] =
                    d_phi_real_squared->item(i, k);
            }
            else
            {
                inverse_input[d_phi_imaginary_squared->numColumns() * 2 * i + j] =
                    d_phi_imaginary_squared->item(i, k);
                k++;
            }
        }
    }

    // Call lapack routines to do the inversion.
    // Set up some stuff the lapack routines need.
    int info;
    int mtx_size = d_phi_real_squared->numColumns();
    int lwork = mtx_size*mtx_size*std::max(10,d_num_procs);
    int* ipiv = new int [mtx_size];
    double* work = new double [lwork];

    // Now call lapack to do the inversion.
    zgetrf(&mtx_size, &mtx_size, inverse_input, &mtx_size, ipiv, &info);
    zgetri(&mtx_size, inverse_input, &mtx_size, ipiv, work, &lwork, &info);

    d_phi_real_squared_inverse = d_phi_real_squared;
    d_phi_imaginary_squared_inverse = d_phi_imaginary_squared;

    for (int i = 0; i < d_phi_real_squared_inverse->numRows(); i++)
    {
        int k = 0;
        for (int j = 0; j < d_phi_real_squared_inverse->numColumns() * 2; j++)
        {
            if (j % 2 == 0)
            {
                d_phi_real_squared_inverse->item(i, k) =
                    inverse_input[d_phi_real_squared_inverse->numColumns() * 2 * i + j];
            }
            else
            {
                d_phi_imaginary_squared_inverse->item(i, k) =
                    inverse_input[d_phi_imaginary_squared_inverse->numColumns() * 2 * i + j];
                k++;
            }
        }
    }

    // Initial condition
    std::unique_ptr<Vector> init_real = d_phi_real->transposeMult(init);
    std::unique_ptr<Vector> init_imaginary = d_phi_imaginary->transposeMult(init);

    std::unique_ptr<Vector> d_projected_init_real_1 =
        d_phi_real_squared_inverse->mult(*init_real);
    std::unique_ptr<Vector> d_projected_init_real_2 =
        d_phi_imaginary_squared_inverse->mult(*init_imaginary);
    d_projected_init_real = d_projected_init_real_1->plus(*d_projected_init_real_2);

    std::unique_ptr<Vector> d_projected_init_imaginary_1 =
        d_phi_real_squared_inverse->mult(*init_imaginary);
    std::unique_ptr<Vector> d_projected_init_imaginary_2 =
        d_phi_imaginary_squared_inverse->mult(*init_real);
    d_projected_init_imaginary = d_projected_init_imaginary_2->minus(
                                     *d_projected_init_imaginary_1);

    // Controls
    std::unique_ptr<Matrix> B_tilde_f = d_B_tilde->mult(controls);
    std::unique_ptr<Matrix> UBf = d_basis->mult(*B_tilde_f);
    std::unique_ptr<Matrix> controls_real = d_phi_real->transposeMult(*UBf);
    std::unique_ptr<Matrix> controls_imaginary = d_phi_imaginary->transposeMult(
                *UBf);

    d_projected_controls_real = d_phi_real_squared_inverse->mult(*controls_real);
    std::unique_ptr<Matrix> d_projected_controls_real_2 =
        d_phi_imaginary_squared_inverse->mult(*controls_imaginary);
    *d_projected_controls_real += *d_projected_controls_real_2;

    d_projected_controls_imaginary = d_phi_imaginary_squared_inverse->mult(
                                         *controls_real);
    std::unique_ptr<Matrix> d_projected_controls_imaginary_2 =
        d_phi_real_squared_inverse->mult(*controls_imaginary);
    *d_projected_controls_imaginary -= *d_projected_controls_imaginary_2;

    delete [] inverse_input;
    delete [] ipiv;
    delete [] work;

    if (t_offset >= 0.0)
    {
        std::cout << "t_offset is updated from " << d_t_offset
                  << " to " << t_offset << std::endl;
        d_t_offset = t_offset;
    }
    d_init_projected = true;
}

std::unique_ptr<Vector>
DMDc::predict(double t)
{
    CAROM_VERIFY(d_trained);
    CAROM_VERIFY(d_init_projected);
    CAROM_VERIFY(t >= 0.0);

    t -= d_t_offset;

    int n = round(t / d_dt);
    //int n = min(round(t / d_dt), d_projected_controls_real->numColumns());

    std::pair<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>> d_phi_pair =
                phiMultEigs(t);
    std::shared_ptr<Matrix> d_phi_mult_eigs_real = d_phi_pair.first;
    std::shared_ptr<Matrix> d_phi_mult_eigs_imaginary = d_phi_pair.second;

    std::unique_ptr<Vector> d_predicted_state_real_1 = d_phi_mult_eigs_real->mult(
                *d_projected_init_real);
    std::unique_ptr<Vector> d_predicted_state_real_2 =
        d_phi_mult_eigs_imaginary->mult(*d_projected_init_imaginary);
    std::unique_ptr<Vector> d_predicted_state_real =
        d_predicted_state_real_1->minus(*d_predicted_state_real_2);
    addOffset(*d_predicted_state_real);

    Vector* f_control_real = new Vector(d_basis->numRows(), false);
    Vector* f_control_imaginary = new Vector(d_basis->numRows(), false);
    for (int k = 0; k < n; k++)
    {
        t -= d_dt;
        std::pair<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>> d_phi_pair =
                    phiMultEigs(t);
        d_phi_mult_eigs_real = d_phi_pair.first;
        d_phi_mult_eigs_imaginary = d_phi_pair.second;

        d_projected_controls_real->getColumn(k, *f_control_real);
        d_projected_controls_imaginary->getColumn(k, *f_control_imaginary);
        d_predicted_state_real_1 = d_phi_mult_eigs_real->mult(
                                       *f_control_real);
        d_predicted_state_real_2 = d_phi_mult_eigs_imaginary->mult(
                                       *f_control_imaginary);
        *d_predicted_state_real += *d_predicted_state_real_1;
        *d_predicted_state_real -= *d_predicted_state_real_2;
    }

    delete f_control_real;
    delete f_control_imaginary;

    return d_predicted_state_real;
}

void
DMDc::addOffset(Vector & result)
{
    if (d_state_offset)
    {
        result += *d_state_offset;
    }
}

std::complex<double>
DMDc::computeEigExp(std::complex<double> eig, double t)
{
    return std::pow(eig, t / d_dt);
}

std::pair<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>>
        DMDc::phiMultEigs(double t)
{
    Matrix* d_eigs_exp_real = new Matrix(d_k, d_k, false);
    Matrix* d_eigs_exp_imaginary = new Matrix(d_k, d_k, false);

    for (int i = 0; i < d_k; i++)
    {
        std::complex<double> eig_exp = computeEigExp(d_eigs[i], t);
        d_eigs_exp_real->item(i, i) = std::real(eig_exp);
        d_eigs_exp_imaginary->item(i, i) = std::imag(eig_exp);
    }

    std::shared_ptr<Matrix> d_phi_mult_eigs_real = d_phi_real->mult(
                *d_eigs_exp_real);
    std::unique_ptr<Matrix> d_phi_mult_eigs_real_2 = d_phi_imaginary->mult(
                *d_eigs_exp_imaginary);
    *d_phi_mult_eigs_real -= *d_phi_mult_eigs_real_2;
    std::shared_ptr<Matrix> d_phi_mult_eigs_imaginary = d_phi_real->mult(
                *d_eigs_exp_imaginary);
    std::unique_ptr<Matrix> d_phi_mult_eigs_imaginary_2 = d_phi_imaginary->mult(
                *d_eigs_exp_real);
    *d_phi_mult_eigs_imaginary += *d_phi_mult_eigs_imaginary_2;

    delete d_eigs_exp_real;
    delete d_eigs_exp_imaginary;

    return std::pair<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>>
           (d_phi_mult_eigs_real,
            d_phi_mult_eigs_imaginary);
}

double
DMDc::getTimeOffset() const
{
    return d_t_offset;
}

std::unique_ptr<const Matrix>
DMDc::getSnapshotMatrix()
{
    return createSnapshotMatrix(d_snapshots);
}

std::unique_ptr<const Matrix>
DMDc::createSnapshotMatrix(const std::vector<std::shared_ptr<Vector>> &
                           snapshots)
{
    CAROM_VERIFY(snapshots.size() > 0);
    CAROM_VERIFY(snapshots[0]->dim() > 0);
    for (int i = 0 ; i < snapshots.size() - 1; i++)
    {
        CAROM_VERIFY(snapshots[i]->dim() == snapshots[i + 1]->dim());
        CAROM_VERIFY(snapshots[i]->distributed() == snapshots[i + 1]->distributed());
    }

    Matrix* snapshot_mat = new Matrix(snapshots[0]->dim(), snapshots.size(),
                                      snapshots[0]->distributed());

    for (int i = 0; i < snapshots[0]->dim(); i++)
    {
        for (int j = 0; j < snapshots.size(); j++)
        {
            snapshot_mat->item(i, j) = snapshots[j]->item(i);
        }
    }

    return std::unique_ptr<const Matrix>(snapshot_mat);
}

void
DMDc::load(std::string base_file_name)
{
    CAROM_ASSERT(!base_file_name.empty());

    char tmp[100];
    std::string full_file_name = base_file_name;
    HDFDatabase database;
    database.open(full_file_name, "r");

    sprintf(tmp, "dt");
    database.getDouble(tmp, d_dt);

    sprintf(tmp, "t_offset");
    database.getDouble(tmp, d_t_offset);

    sprintf(tmp, "k");
    database.getInteger(tmp, d_k);

    sprintf(tmp, "num_eigs");
    int num_eigs;
    database.getInteger(tmp, num_eigs);

    std::vector<double> eigs_real;
    std::vector<double> eigs_imag;

    sprintf(tmp, "eigs_real");
    eigs_real.resize(num_eigs);
    database.getDoubleArray(tmp, &eigs_real[0], num_eigs);

    sprintf(tmp, "eigs_imag");
    eigs_imag.resize(num_eigs);
    database.getDoubleArray(tmp, &eigs_imag[0], num_eigs);

    for (int i = 0; i < num_eigs; i++)
    {
        d_eigs.push_back(std::complex<double>(eigs_real[i], eigs_imag[i]));
    }
    database.close();

    full_file_name = base_file_name + "_basis";
    d_basis.reset(new Matrix());
    d_basis->read(full_file_name);

    full_file_name = base_file_name + "_A_tilde";
    d_A_tilde.reset(new Matrix());
    d_A_tilde->read(full_file_name);

    full_file_name = base_file_name + "_B_tilde";
    d_B_tilde.reset(new Matrix());
    d_B_tilde->read(full_file_name);

    full_file_name = base_file_name + "_phi_real";
    d_phi_real.reset(new Matrix());
    d_phi_real->read(full_file_name);

    full_file_name = base_file_name + "_phi_imaginary";
    d_phi_imaginary.reset(new Matrix());
    d_phi_imaginary->read(full_file_name);

    full_file_name = base_file_name + "_phi_real_squared_inverse";
    d_phi_real_squared_inverse.reset(new Matrix());
    d_phi_real_squared_inverse->read(full_file_name);

    full_file_name = base_file_name + "_phi_imaginary_squared_inverse";
    d_phi_imaginary_squared_inverse.reset(new Matrix());
    d_phi_imaginary_squared_inverse->read(full_file_name);

    full_file_name = base_file_name + "_projected_init_real";
    d_projected_init_real.reset(new Vector());
    d_projected_init_real->read(full_file_name);

    full_file_name = base_file_name + "_projected_init_imaginary";
    d_projected_init_imaginary.reset(new Vector());
    d_projected_init_imaginary->read(full_file_name);

    full_file_name = base_file_name + "_state_offset";
    if (Utilities::file_exist(full_file_name + ".000000"))
    {
        d_state_offset.reset(new Vector());
        d_state_offset->read(full_file_name);
    }

    d_init_projected = true;
    d_trained = true;

    MPI_Barrier(MPI_COMM_WORLD);
}

void
DMDc::load(const char* base_file_name)
{
    load(std::string(base_file_name));
}

void
DMDc::save(std::string base_file_name)
{
    CAROM_ASSERT(!base_file_name.empty());
    CAROM_VERIFY(d_trained);

    if (d_rank == 0)
    {
        char tmp[100];
        std::string full_file_name = base_file_name;
        HDFDatabase database;
        database.create(full_file_name);

        sprintf(tmp, "dt");
        database.putDouble(tmp, d_dt);

        sprintf(tmp, "t_offset");
        database.putDouble(tmp, d_t_offset);

        sprintf(tmp, "k");
        database.putInteger(tmp, d_k);

        sprintf(tmp, "num_eigs");
        database.putInteger(tmp, d_eigs.size());

        std::vector<double> eigs_real;
        std::vector<double> eigs_imag;

        for (int i = 0; i < d_eigs.size(); i++)
        {
            eigs_real.push_back(d_eigs[i].real());
            eigs_imag.push_back(d_eigs[i].imag());
        }

        sprintf(tmp, "eigs_real");
        database.putDoubleArray(tmp, &eigs_real[0], eigs_real.size());

        sprintf(tmp, "eigs_imag");
        database.putDoubleArray(tmp, &eigs_imag[0], eigs_imag.size());
        database.close();
    }

    std::string full_file_name;

    if (d_basis != NULL)
    {
        full_file_name = base_file_name + "_basis";
        d_basis->write(full_file_name);
    }

    if (d_A_tilde != NULL)
    {
        full_file_name = base_file_name + "_A_tilde";
        d_A_tilde->write(full_file_name);
    }

    if (d_B_tilde != NULL)
    {
        full_file_name = base_file_name + "_B_tilde";
        d_B_tilde->write(full_file_name);
    }

    full_file_name = base_file_name + "_phi_real";
    d_phi_real->write(full_file_name);

    full_file_name = base_file_name + "_phi_imaginary";
    d_phi_imaginary->write(full_file_name);

    full_file_name = base_file_name + "_phi_real_squared_inverse";
    d_phi_real_squared_inverse->write(full_file_name);

    full_file_name = base_file_name + "_phi_imaginary_squared_inverse";
    d_phi_imaginary_squared_inverse->write(full_file_name);

    full_file_name = base_file_name + "_projected_init_real";
    d_projected_init_real->write(full_file_name);

    full_file_name = base_file_name + "_projected_init_imaginary";
    d_projected_init_imaginary->write(full_file_name);

    if (d_state_offset != NULL)
    {
        full_file_name = base_file_name + "_state_offset";
        d_state_offset->write(full_file_name);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

void
DMDc::save(const char* base_file_name)
{
    save(std::string(base_file_name));
}

void
DMDc::summary(std::string base_file_name)
{
    if (d_rank == 0)
    {
        CSVDatabase* csv_db(new CSVDatabase);

        csv_db->putDoubleVector(base_file_name + "_singular_value.csv", d_sv,
                                d_num_singular_vectors);
        csv_db->putComplexVector(base_file_name + "_eigenvalue.csv", d_eigs,
                                 d_eigs.size());

        delete csv_db;
    }
}

}
