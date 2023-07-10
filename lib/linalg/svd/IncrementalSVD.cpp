/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The abstract incremental SVD algorithm defines algorithm
//              interface.

#include "IncrementalSVD.h"
#include "utils/HDFDatabase.h"

#include "mpi.h"

#include <cmath>
#include <iomanip>
#include <limits>
#include <stdio.h>
#include <sstream>

/* Use Autotools-detected Fortran name-mangling scheme */
#define dgesdd CAROM_FC_GLOBAL(dgesdd, DGESDD)

extern "C" {
    void dgesdd(char*, int*, int*, double*, int*,
                double*, double*, int*, double*, int*,
                double*, int*, int*, int*);
}

namespace CAROM {

const int IncrementalSVD::COMMUNICATE_U = 666;

IncrementalSVD::IncrementalSVD(
    Options options,
    const std::string& basis_file_name) :
    SVD(options),
    d_linearity_tol(options.linearity_tol),
    d_skip_linearly_dependent(options.skip_linearly_dependent),
    d_max_basis_dimension(options.max_basis_dimension),
    d_total_dim(0),
    d_save_state(options.save_state),
    d_update_right_SV(options.update_right_SV),
    d_state_database(0)
{
    CAROM_VERIFY(options.linearity_tol > 0.0);
    CAROM_VERIFY(options.max_basis_dimension > 0);

    // Get the number of processors, the dimensions for each process, and the
    // total dimension.
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init) {
        MPI_Comm_size(MPI_COMM_WORLD, &d_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    }
    else {
        d_size = 1;
        d_rank = 0;
    }
    d_proc_dims.reserve(d_size);
    d_proc_dims.resize(d_size);
    if (mpi_init) {
        MPI_Allgather(&d_dim,
                      1,
                      MPI_INT,
                      &d_proc_dims[0],
                      1,
                      MPI_INT,
                      MPI_COMM_WORLD);
    }
    else {
        d_proc_dims[0] = d_dim;
    }
    for (int i = 0; i < d_size; ++i) {
        d_total_dim += d_proc_dims[i];
    }

    // If the state of the SVD is to be restored then open the database and
    // restore the necessary data from the database now.
    if (options.save_state || options.restore_state) {
        std::ostringstream tmp;
        tmp << basis_file_name << ".state." <<
            std::setw(6) << std::setfill('0') << d_rank;
        d_state_file_name = tmp.str();
    }
    if (options.restore_state) {
        // Open state database file.
        d_state_database = new HDFDatabase();
        bool is_good = d_state_database->open(d_state_file_name, "r");
        if (is_good) {
            // Read time interval start time.
            double time;
            d_state_database->getDouble("time", time);
            d_time_interval_start_times.resize(1);
            d_time_interval_start_times[0] = time;

            // Read d_U.
            int num_rows;
            d_state_database->getInteger("U_num_rows", num_rows);
            int num_cols;
            d_state_database->getInteger("U_num_cols", num_cols);
            d_U = new Matrix(num_rows, num_cols, true);
            d_state_database->getDoubleArray("U",
                                             &d_U->item(0, 0),
                                             num_rows*num_cols);

            if (d_update_right_SV) {
                // Read d_W.
                d_state_database->getInteger("W_num_rows", num_rows);
                d_state_database->getInteger("W_num_cols", num_cols);
                d_W = new Matrix(num_rows, num_cols, true);
                d_state_database->getDoubleArray("W",
                                                 &d_W->item(0, 0),
                                                 num_rows*num_cols);
            }

            // Read d_S.
            int num_dim;
            d_state_database->getInteger("S_num_dim", num_dim);
            d_S = new Vector(num_dim, false);
            d_state_database->getDoubleArray("S",
                                             &d_S->item(0),
                                             num_dim);

            // Set d_num_samples.
            d_num_samples = num_cols;
        }
        else {
            delete d_state_database;
            d_state_database = 0;
        }
    }
}

IncrementalSVD::~IncrementalSVD()
{
    // If the state of the SVD is to be saved, then save d_S and d_U now.  The
    // derived class has already created the database.
    //
    // If there are multiple time intervals then saving and restoring the state
    // does not make sense as there is not one, all encompassing, basis.
    if (d_save_state && d_time_interval_start_times.size() == 1) {
        // Save the time interval start time.
        d_state_database->putDouble("time", d_time_interval_start_times[0]);

        // Save d_U.
        int num_rows = d_U->numRows();
        d_state_database->putInteger("U_num_rows", num_rows);
        int num_cols = d_U->numColumns();
        d_state_database->putInteger("U_num_cols", num_cols);
        d_state_database->putDoubleArray("U", &d_U->item(0, 0), num_rows*num_cols);

        // Save d_S.
        int num_dim = d_S->dim();
        d_state_database->putInteger("S_num_dim", num_dim);
        d_state_database->putDoubleArray("S",
                                         &d_S->item(0),
                                         num_dim);

        if (d_update_right_SV) {
            // Save d_W.
            num_rows = d_W->numRows();
            d_state_database->putInteger("W_num_rows", num_rows);
            num_cols = d_W->numColumns();
            d_state_database->putInteger("W_num_cols", num_cols);
            d_state_database->putDoubleArray("W",
                                             &d_W->item(0, 0),
                                             num_rows*num_cols);
        }

        // Close state database file and delete database object.
        d_state_database->close();
        delete d_state_database;
    }
}

bool
IncrementalSVD::takeSample(
    double* u_in,
    double time,
    bool add_without_increase)
{
    CAROM_VERIFY(u_in != 0);
    CAROM_VERIFY(time >= 0.0);

    // Check that u_in is not non-zero.
    Vector u_vec(u_in, d_dim, true);
    if (u_vec.norm() == 0.0) {
        return false;
    }

    // If this is the first SVD then build it.  Otherwise add this sample to the
    // system.
    bool result = true;
    if (isNewTimeInterval()) {
        buildInitialSVD(u_in, time);
    }
    else {
        result = buildIncrementalSVD(u_in,add_without_increase);
    }

    if (d_debug_algorithm) {
        const Matrix* basis = getSpatialBasis();
        if (d_rank == 0) {
            // Print d_S.
            for (int col = 0; col < d_num_samples; ++col) {
                printf("%.16e  ", d_S->item(col));
                printf("\n");
            }
            printf("\n");

            // Print process 0's part of the basis.
            for (int row = 0; row < d_dim; ++row) {
                for (int col = 0; col < d_num_samples; ++col) {
                    printf("%.16e ", basis->item(row, col));
                }
                printf("\n");
            }

            // Gather other processor's parts of the basis and print them.
            for (int proc = 1; proc < d_size; ++proc) {
                double* m = new double[d_proc_dims[proc]*d_num_samples];
                MPI_Status status;
                MPI_Recv(m,
                         d_proc_dims[proc]*d_num_samples,
                         MPI_DOUBLE,
                         proc,
                         COMMUNICATE_U,
                         MPI_COMM_WORLD,
                         &status);
                int idx = 0;
                for (int row = 0; row < d_proc_dims[proc]; ++row) {
                    for (int col = 0; col < d_num_samples; ++col) {
                        printf("%.16e ", m[idx++]);
                    }
                    printf("\n");
                }
                delete [] m;
            }
            printf("============================================================\n");
        }
        else {
            // Send this processor's part of the basis to process 0.
            MPI_Request request;
            MPI_Isend(const_cast<double*>(&basis->item(0, 0)),
                      d_dim*d_num_samples,
                      MPI_DOUBLE,
                      0,
                      COMMUNICATE_U,
                      MPI_COMM_WORLD,
                      &request);
        }
    }
    return result;
}

const Matrix*
IncrementalSVD::getSpatialBasis()
{
    CAROM_ASSERT(d_basis != 0);
    return d_basis;
}

const Matrix*
IncrementalSVD::getTemporalBasis()
{
    CAROM_ASSERT(d_basis != 0);
    return d_basis_right;
}

const Vector*
IncrementalSVD::getSingularValues()
{
    CAROM_ASSERT(d_S != 0);
    return d_S;
}

const Matrix*
IncrementalSVD::getSnapshotMatrix()
{
    std::cout << "Getting snapshot matrix for incremental SVD not implemented" <<
              std::endl;
    return d_snapshots;
}

bool
IncrementalSVD::buildIncrementalSVD(
    double* u, bool add_without_increase)
{
    CAROM_VERIFY(u != 0);

    // l = basis' * u
    Vector u_vec(u, d_dim, true);
    Vector* l = d_basis->transposeMult(u_vec);

    // basisl = basis * l
    Vector* basisl = d_basis->mult(l);

    // Computing as k = sqrt(u.u - 2.0*l.l + basisl.basisl)
    // results in catastrophic cancellation, and must be avoided.
    // Instead we compute as k = sqrt((u-basisl).(u-basisl)).
    Vector* e_proj = u_vec.minus(basisl);
    double k = e_proj->inner_product(e_proj);
    delete e_proj;

    if (k <= 0) {
        if(d_rank == 0) printf("linearly dependent sample!\n");
        k = 0;
    }
    else {
        k = sqrt(k);
    }

    // Use k to see if the vector addressed by u is linearly dependent
    // on the left singular vectors.
    bool linearly_dependent_sample;
    if ( k < d_linearity_tol ) {
        if(d_rank == 0) {
            std::cout << "linearly dependent sample! k = " << k << "\n";
            std::cout << "d_linearity_tol = " << d_linearity_tol << "\n";
        }
        k = 0;
        linearly_dependent_sample = true;
    } else if ( d_num_samples >= d_max_basis_dimension || add_without_increase ) {
        k = 0;
        linearly_dependent_sample = true;
    }
    // Check to see if the "number of samples" (in IncrementalSVD and
    // its subclasses, d_num_samples appears to be equal to the number
    // of columns of the left singular vectors) is greater than or equal
    // to the dimension of snapshot vectors. If so, then the vector
    // addressed by the pointer u must be linearly dependent on the left
    // singular vectors.
    else if (d_num_samples >= d_total_dim) {
        linearly_dependent_sample = true;
    }
    else {
        linearly_dependent_sample = false;
    }

    // Create Q.
    double* Q;
    constructQ(Q, l, k);
    delete l;

    // Now get the singular value decomposition of Q.
    Matrix* A;
    Matrix* W;
    Matrix* sigma;
    bool result = svd(Q, A, sigma, W);

    // Done with Q.
    delete [] Q;

    // If the svd was successful then add the sample.  Otherwise clean up and
    // return.
    if (result) {

        // We need to add the sample if it is not linearly dependent or if it is
        // linearly dependent and we are not skipping linearly dependent samples.
        if ((linearly_dependent_sample && !d_skip_linearly_dependent) ) {
            // This sample is linearly dependent and we are not skipping linearly
            // dependent samples.
            if(d_rank == 0) std::cout << "adding linearly dependent sample!\n";
            addLinearlyDependentSample(A, W, sigma);
            delete sigma;
        }
        else if (!linearly_dependent_sample) {
            // This sample is not linearly dependent.

            // Compute j
            Vector* j = u_vec.minus(basisl);
            for (int i = 0; i < d_dim; ++i) {
                j->item(i) /= k;
            }

            // addNewSample will assign sigma to d_S hence it should not be
            // deleted upon return.
            addNewSample(j, A, W, sigma);
            delete j;
        }
        delete basisl;
        delete A;
        delete W;

        // Compute the basis vectors.
        computeBasis();
    }
    else {
        delete basisl;
        delete A;
        delete W;
        delete sigma;
    }
    return result;
}

void
IncrementalSVD::constructQ(
    double*& Q,
    const Vector* l,
    double k)
{
    CAROM_VERIFY(l != 0);
    CAROM_VERIFY(l->dim() == numSamples());

    // Create Q.
    Q = new double [(d_num_samples+1)*(d_num_samples+1)];

    // Fill Q in column major order.
    int q_idx = 0;
    for (int row = 0; row < d_num_samples; ++row) {
        q_idx = row;
        for (int col = 0; col < d_num_samples; ++col) {
            if (row == col)
            {
                Q[q_idx] = d_S->item(col);
            }
            else
            {
                Q[q_idx] = 0.0;
            }
            q_idx += d_num_samples+1;
        }
        Q[q_idx] = l->item(row);
    }
    q_idx = d_num_samples;
    for (int col = 0; col < d_num_samples; ++col) {
        Q[q_idx] = 0.0;
        q_idx += d_num_samples+1;
    }
    Q[q_idx] = k;
}

bool
IncrementalSVD::svd(
    double* A,
    Matrix*& U,
    Matrix*& S,
    Matrix*& V)
{
    CAROM_VERIFY(A != 0);

    // Construct U, S, and V.
    U = new Matrix(d_num_samples+1, d_num_samples+1, false);
    S = new Matrix(d_num_samples+1, d_num_samples+1, false);
    V = new Matrix(d_num_samples+1, d_num_samples+1, false);
    for (int row = 0; row < d_num_samples+1; ++row) {
        for (int col = 0; col < d_num_samples+1; ++col) {
            S->item(row, col) = 0.0;
        }
    }

    // Use lapack's dgesdd Fortran function to perform the svd.  As this is
    // Fortran A and all the computed matrices are in column major order.
    double* sigma = new double [d_num_samples+1];
    char jobz = 'A';
    int m = d_num_samples+1;
    int n = d_num_samples+1;
    int lda = d_num_samples+1;
    int ldu = d_num_samples+1;
    int ldv = d_num_samples+1;
    int lwork = m*(4*m + 7);
    double* work = new double [lwork];
    int iwork[8*m];
    int info;
    dgesdd(&jobz,
           &m,
           &n,
           A,
           &lda,
           sigma,
           &U->item(0, 0),
           &ldu,
           &V->item(0, 0),
           &ldv,
           work,
           &lwork,
           iwork,
           &info);
    delete [] work;

    // If the svd succeeded, fill U and S.  Otherwise clean up and return.
    if (info == 0) {
        // Place sigma into S.
        for (int i = 0; i < d_num_samples+1; ++i) {
            S->item(i, i) = sigma[i];
        }
        delete [] sigma;

        // U is column major order so convert it to row major order.
        for (int row = 0; row < d_num_samples+1; ++row) {
            for (int col = row+1; col < d_num_samples+1; ++col) {
                double tmp = U->item(row, col);
                U->item(row, col) = U->item(col, row);
                U->item(col, row) = tmp;
            }
        }
        /*      if(d_update_right_SV) {
                // V is column major order so convert it to row major order.
                for (int row = 0; row < d_num_samples+1; ++row) {
                   for (int col = row+1; col < d_num_samples+1; ++col) {
                      double tmp = V->item(row, col);
                      V->item(row, col) = V->item(col, row);
                      V->item(col, row) = tmp;
                   }
                }
              }*/
    }
    else {
        delete [] sigma;
    }
    return info == 0;
}

double
IncrementalSVD::checkOrthogonality(
    const Matrix* m)
{
    CAROM_ASSERT(m != 0);

    double result = 0.0;
    if (d_num_samples > 1) {
        int last_col = d_num_samples-1;
        double tmp = 0.0;
        int num_rows = m->numRows();
        for (int i = 0; i < num_rows; ++i) {
            tmp += m->item(i, 0) * m->item(i, last_col);
        }
        if (m->distributed() && d_size > 1) {
            MPI_Allreduce(&tmp, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
        else {
            result = tmp;
        }
    }
    return result;
}

}
