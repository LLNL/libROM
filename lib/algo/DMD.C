/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Implementation of the DMD algorithm.

#include "DMD.h"

#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include "linalg/scalapack_wrapper.h"
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

DMD::DMD(int dim, double dt)
{
    CAROM_VERIFY(dim > 0);
    CAROM_VERIFY(dt > 0);

    // Get the rank of this process, and the number of processors.
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);
    d_dim = dim;
    d_dt = dt;
    d_trained = false;
}

DMD::DMD(std::string base_file_name)
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

    load(base_file_name);
}

void DMD::takeSample(double* u_in, double t)
{
    CAROM_VERIFY(u_in != 0);
    CAROM_VERIFY(t >= 0.0);
    Vector* sample = new Vector(u_in, d_dim, true);
    if (d_snapshots.empty())
    {
        d_t_offset = t;
    }
    d_snapshots.push_back(sample);
}

void DMD::train(double energy_fraction)
{
    const Matrix* f_snapshots = getSnapshotMatrix();
    CAROM_VERIFY(f_snapshots->numColumns() > 1);
    CAROM_VERIFY(energy_fraction > 0 && energy_fraction <= 1);
    d_energy_fraction = energy_fraction;
    constructDMD(f_snapshots, d_rank, d_num_procs);

    delete f_snapshots;
}

void DMD::train(int k)
{
    const Matrix* f_snapshots = getSnapshotMatrix();
    CAROM_VERIFY(f_snapshots->numColumns() > 1);
    CAROM_VERIFY(k > 0 && k <= f_snapshots->numColumns() - 1);
    d_energy_fraction = -1.0;
    d_k = k;
    constructDMD(f_snapshots, d_rank, d_num_procs);

    delete f_snapshots;
}

void
DMD::constructDMD(const Matrix* f_snapshots,
                  int d_rank,
                  int d_num_procs)
{
    // TODO: Making two copies of the snapshot matrix has a lot of overhead.
    //       We need to figure out a way to do submatrix multiplication and to
    //       reimplement this algorithm using one snapshot matrix.
    Matrix* f_snapshots_minus = new Matrix(f_snapshots->numRows(),
                                           f_snapshots->numColumns() - 1, f_snapshots->distributed());
    Matrix* f_snapshots_plus = new Matrix(f_snapshots->numRows(),
                                          f_snapshots->numColumns() - 1, f_snapshots->distributed());

    // Break up snapshots into snapshots_minus and snapshots_plus
    // snapshots_minus = all columns of snapshots except last
    // snapshots_plus = all columns of snapshots except first
    for (int i = 0; i < f_snapshots->numRows(); i++)
    {
        for (int j = 0; j < f_snapshots->numColumns() - 1; j++)
        {
            f_snapshots_minus->item(i, j) = f_snapshots->item(i, j);
            f_snapshots_plus->item(i, j) = f_snapshots->item(i, j + 1);
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
    d_num_singular_vectors = std::min(f_snapshots_minus->numColumns(), f_snapshots_minus->numDistributedRows());
    if (d_energy_fraction != -1.0)
    {
        d_k = d_num_singular_vectors;
        double total_energy = 0.0;
        for (int i = 0; i < d_num_singular_vectors; i++)
        {
            d_sv.push_back(d_factorizer->S[i]);
            total_energy += d_factorizer->S[i];
        }
        double current_energy = 0.0;
        for (int i = 0; i < d_num_singular_vectors; i++)
        {
            current_energy += d_factorizer->S[i];
            if (current_energy / total_energy >= d_energy_fraction)
            {
                d_k = i + 1;
                break;
            }
        }
    }

    std::cout << "Using " << d_k << " basis vectors out of " << d_num_singular_vectors << "." << std::endl;

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

    // Calculate the right eigenvalues/eigenvectors of A_tilde
    ComplexEigenPair eigenpair = NonSymmetricRightEigenSolve(A_tilde);
    d_eigs = eigenpair.eigs;

    // Calculate phi
    Matrix* f_snapshots_plus_mult_d_basis_right = f_snapshots_plus->mult(d_basis_right);
    Matrix* f_snapshots_plus_mult_d_basis_right_mult_d_S_inv = f_snapshots_plus_mult_d_basis_right->mult(d_S_inv);
    d_phi_real = f_snapshots_plus_mult_d_basis_right_mult_d_S_inv->mult(eigenpair.ev_real);
    d_phi_imaginary = f_snapshots_plus_mult_d_basis_right_mult_d_S_inv->mult(eigenpair.ev_imaginary);

    Vector* init = new Vector(f_snapshots_minus->numRows(), true);
    for (int i = 0; i < init->dim(); i++)
    {
        init->item(i) = f_snapshots_minus->item(i, 0);
    }

    // Calculate pinv(d_phi) * initial_condition.
    projectInitialCondition(init);

    d_trained = true;

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
    delete eigenpair.ev_real;
    delete eigenpair.ev_imaginary;
    delete init;
}


void
DMD::projectInitialCondition(const Vector* init)
{
    Matrix* d_phi_real_squared = d_phi_real->transposeMult(d_phi_real);
    Matrix* d_phi_real_squared_2 = d_phi_imaginary->transposeMult(d_phi_imaginary);
    *d_phi_real_squared += *d_phi_real_squared_2;

    Matrix* d_phi_imaginary_squared = d_phi_real->transposeMult(d_phi_imaginary);
    Matrix* d_phi_imaginary_squared_2 = d_phi_imaginary->transposeMult(d_phi_real);
    *d_phi_imaginary_squared -= *d_phi_imaginary_squared_2;

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
    int lwork = mtx_size*mtx_size*std::max(10,d_num_procs);
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

    Vector* d_projected_init_real_1 = d_phi_real_squared->mult(rhs_real);
    Vector* d_projected_init_real_2 = d_phi_imaginary_squared->mult(rhs_imaginary);
    d_projected_init_real = d_projected_init_real_1->plus(d_projected_init_real_2);

    Vector* d_projected_init_imaginary_1 = d_phi_real_squared->mult(rhs_imaginary);
    Vector* d_projected_init_imaginary_2 = d_phi_imaginary_squared->mult(rhs_real);
    d_projected_init_imaginary = d_projected_init_imaginary_2->minus(d_projected_init_imaginary_1);

    delete d_phi_real_squared;
    delete d_phi_real_squared_2;
    delete d_projected_init_real_1;
    delete d_projected_init_real_2;
    delete d_phi_imaginary_squared;
    delete d_phi_imaginary_squared_2;
    delete d_projected_init_imaginary_1;
    delete d_projected_init_imaginary_2;
    delete rhs_real;
    delete rhs_imaginary;

    delete [] inverse_input;
    delete [] ipiv;
    delete [] work;
}

Vector*
DMD::predict(double t)
{
    const std::pair<Vector*, Vector*> d_projected_init_pair(d_projected_init_real, d_projected_init_imaginary);
    return predict(d_projected_init_pair, t);
}

Vector*
DMD::predict(const std::pair<Vector*, Vector*> init,
             double t)
{
    CAROM_VERIFY(d_trained);
    CAROM_VERIFY(t >= 0.0);
    t -= d_t_offset;
    double t_dt = t / d_dt;
    std::pair<Matrix*, Matrix*> d_phi_pair = phiMultEigs(t_dt);
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
DMD::phiMultEigs(double t)
{
    Matrix* d_eigs_exp_real = new Matrix(d_k, d_k, false);
    Matrix* d_eigs_exp_imaginary = new Matrix(d_k, d_k, false);

    for (int i = 0; i < d_k; i++)
    {
        std::complex<double> eig_exp = std::pow(d_eigs[i], t);
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

double
DMD::getTimeOffset() const
{
    return d_t_offset;
}

const Matrix*
DMD::getSnapshotMatrix()
{
    return createSnapshotMatrix(d_snapshots);
}

const Matrix*
DMD::createSnapshotMatrix(std::vector<Vector*> snapshots)
{
    CAROM_VERIFY(snapshots.size() > 0);
    CAROM_VERIFY(snapshots[0]->dim() > 0);
    for (int i = 0 ; i < snapshots.size() - 1; i++)
    {
        CAROM_VERIFY(snapshots[i]->dim() == snapshots[i + 1]->dim());
        CAROM_VERIFY(snapshots[i]->distributed() == snapshots[i + 1]->distributed());
    }

    Matrix* snapshot_mat = new Matrix(snapshots[0]->dim(), snapshots.size(), snapshots[0]->distributed());

    for (int i = 0; i < snapshots[0]->dim(); i++)
    {
        for (int j = 0; j < snapshots.size(); j++)
        {
            snapshot_mat->item(i, j) = snapshots[j]->item(i);
        }
    }

    return snapshot_mat;
}

void
DMD::load(std::string base_file_name)
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

    full_file_name = base_file_name + "_phi_real";
    d_phi_real = new Matrix();
    d_phi_real->read(full_file_name);

    full_file_name = base_file_name + "_phi_imaginary";
    d_phi_imaginary = new Matrix();
    d_phi_imaginary->read(full_file_name);

    full_file_name = base_file_name + "_projected_init_real";
    d_projected_init_real = new Vector();
    d_projected_init_real->read(full_file_name);

    full_file_name = base_file_name + "_projected_init_imaginary";
    d_projected_init_imaginary = new Vector();
    d_projected_init_imaginary->read(full_file_name);

    MPI_Barrier(MPI_COMM_WORLD);
}

void
DMD::save(std::string base_file_name)
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

    std::string full_file_name = base_file_name + "_phi_real";
    d_phi_real->write(full_file_name);

    full_file_name = base_file_name + "_phi_imaginary";
    d_phi_imaginary->write(full_file_name);

    full_file_name = base_file_name + "_projected_init_real";
    d_projected_init_real->write(full_file_name);

    full_file_name = base_file_name + "_projected_init_imaginary";
    d_projected_init_imaginary->write(full_file_name);

    MPI_Barrier(MPI_COMM_WORLD);
}

void
DMD::summary(std::string output_path)
{
    if (d_rank == 0)
    {
        CSVDatabase* csv_db(new CSVDatabase);

        csv_db->putDoubleVector(output_path + "/singular_value.csv", d_sv, d_num_singular_vectors);
        csv_db->putComplexVector(output_path + "/eigenvalue.csv", d_eigs, d_eigs.size());

        delete csv_db;
    }
}

}
