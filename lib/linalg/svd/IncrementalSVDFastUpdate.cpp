/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The concrete implementation of the incremental SVD algorithm
//              using Matthew Brand's "fast update" method.

#include "IncrementalSVDFastUpdate.h"
#include "utils/HDFDatabase.h"

#include "mpi.h"

#include <cmath>
#include <limits>

namespace CAROM {

IncrementalSVDFastUpdate::IncrementalSVDFastUpdate(
    Options options,
    const std::string& basis_file_name) :
    IncrementalSVD(
        options,
        basis_file_name),
    d_Up(0),
    d_singular_value_tol(options.singular_value_tol)
{
    CAROM_VERIFY(options.singular_value_tol >= 0);

    // If the state of the SVD is to be restored, do it now.  The base class,
    // IncrementalSVD, has already opened the database and restored the state
    // common to all incremental algorithms.  This particular class must also
    // read the state of d_Up and then compute the basis.  If the database could
    // not be found then we can not restore the state.
    if (options.restore_state && d_state_database) {
        // Read d_Up.
        int num_rows;
        d_state_database->getInteger("Up_num_rows", num_rows);
        int num_cols;
        d_state_database->getInteger("Up_num_cols", num_cols);
        d_Up = new Matrix(num_rows, num_cols, true);
        d_state_database->getDoubleArray("Up",
                                         &d_Up->item(0, 0),
                                         num_rows*num_cols);

        // Close and delete the database.
        d_state_database->close();
        delete d_state_database;

        // Compute the basis.
        computeBasis();
    }
}

IncrementalSVDFastUpdate::~IncrementalSVDFastUpdate()
{
    // If the state of the SVD is to be saved, then create the database now.
    // The IncrementalSVD base class destructor will save d_S and d_U.  This
    // derived class must save its specific state data, d_Up.
    //
    // If there are multiple time intervals then saving and restoring the state
    // does not make sense as there is not one, all encompassing, basis.
    if (d_save_state && d_time_interval_start_times.size() == 1) {
        // Create state database file.
        d_state_database = new HDFDatabase();
        d_state_database->create(d_state_file_name);

        // Save d_Up.
        int num_rows = d_Up->numRows();
        d_state_database->putInteger("Up_num_rows", num_rows);
        int num_cols = d_Up->numColumns();
        d_state_database->putInteger("Up_num_cols", num_cols);
        d_state_database->putDoubleArray("Up",
                                         &d_Up->item(0, 0),
                                         num_rows*num_cols);
    }

    // Delete data members.
    if (d_Up) {
        delete d_Up;
    }
}

void
IncrementalSVDFastUpdate::buildInitialSVD(
    double* u,
    double time)
{
    CAROM_VERIFY(u != 0);
    CAROM_VERIFY(time >= 0.0);

    // We have a new time interval.

    // If this is not the first time interval then delete d_basis, d_U, d_Up,
    // and d_S of the just completed time interval.
    int num_time_intervals =
        static_cast<int>(d_time_interval_start_times.size());
    if (num_time_intervals > 0) {
        delete d_basis;
        delete d_U;
        delete d_Up;
        delete d_S;
        delete d_W;
    }
    increaseTimeInterval();
    d_time_interval_start_times[num_time_intervals] = time;

    // Build d_S for this new time interval.
    d_S = new Vector(1, false);
    Vector u_vec(u, d_dim, true);
    double norm_u = u_vec.norm();
    d_S->item(0) = norm_u;

    // Build d_Up for this new time interval.
    d_Up = new Matrix(1, 1, false);
    d_Up->item(0, 0) = 1.0;

    // Build d_U for this new time interval.
    d_U = new Matrix(d_dim, 1, true);
    for (int i = 0; i < d_dim; ++i) {
        d_U->item(i, 0) = u[i]/norm_u;
    }

    // Build d_W for this new time interval.
    if (d_update_right_SV) {
        d_W = new Matrix(1, 1, false);
        d_W->item(0, 0) = 1.0;
    }

    // We now have the first sample for the new time interval.
    d_num_samples = 1;
    d_num_rows_of_W = 1;

    // Compute the basis vectors for this time interval.
    computeBasis();

}

void
IncrementalSVDFastUpdate::computeBasis()
{
    d_basis = d_U->mult(d_Up);
    if(d_update_right_SV)
    {
        delete d_basis_right;
        d_basis_right = new Matrix(*d_W);
    }

    if(d_rank == 0) {
        std::cout << "d_num_samples = " << d_num_samples << "\n";
        std::cout << "d_num_rows_of_W = " << d_num_rows_of_W << "\n";
        std::cout << "d_singular_value_tol = " << d_singular_value_tol << "\n";
        std::cout << "smallest SV = " << d_S->item(d_num_samples-1) << "\n";
        if (d_num_samples > 1) {
            std::cout << "next smallest SV = " << d_S->item(d_num_samples-2) << "\n";
        }
    }
    // remove the smallest singular value if it is smaller than d_singular_value_tol
    if ( (d_singular_value_tol != 0.0) &&
            (d_S->item(d_num_samples-1) < d_singular_value_tol) &&
            (d_num_samples != 1) ) {

        if (d_rank == 0) std::cout << "removing a small singular value!\n";

        Matrix* d_basis_new = new Matrix(d_dim, d_num_samples-1,
                                         d_basis->distributed());
        for (int row = 0; row < d_dim; ++row) {
            for (int col = 0; col < d_num_samples-1; ++col) {
                d_basis_new->item(row, col) = d_basis->item(row,col);
            }
        }
        delete d_basis;
        d_basis = d_basis_new;

        if (d_update_right_SV)
        {
            Matrix* d_basis_right_new = new Matrix(d_num_rows_of_W, d_num_samples-1,
                                                   d_basis_right->distributed());
            for (int row = 0; row < d_num_rows_of_W; ++row) {
                for (int col = 0; col < d_num_samples-1; ++col) {
                    d_basis_right_new->item(row, col) = d_basis_right->item(row,col);
                }
            }
            delete d_basis_right;
            d_basis_right = d_basis_right_new;
        }
        --d_num_samples;
    }

    // Reorthogonalize if necessary.
    if (fabs(checkOrthogonality(d_basis)) >
            std::numeric_limits<double>::epsilon()*static_cast<double>(d_num_samples)) {
        d_basis->orthogonalize();
    }
    if(d_update_right_SV)
    {
        if (fabs(checkOrthogonality(d_basis_right)) >
                std::numeric_limits<double>::epsilon()*d_num_samples) {
            d_basis_right->orthogonalize();
        }
    }

}

void
IncrementalSVDFastUpdate::addLinearlyDependentSample(
    const Matrix* A,
    const Matrix* W,
    const Matrix* sigma)
{
    CAROM_VERIFY(A != 0);
    CAROM_VERIFY(sigma != 0);

    // Chop a row and a column off of A to form Amod.  Also form
    // d_S by chopping a row and a column off of sigma.
    Matrix Amod(d_num_samples, d_num_samples, false);
    for (int row = 0; row < d_num_samples; ++row) {
        for (int col = 0; col < d_num_samples; ++col) {
            Amod.item(row, col) = A->item(row, col);
            if (row == col)
            {
                d_S->item(col) = sigma->item(row, col);
            }
        }
    }

    // Multiply d_Up and Amod and put result into d_Up.
    Matrix* Up_times_Amod = d_Up->mult(Amod);
    delete d_Up;
    d_Up = Up_times_Amod;

    Matrix* new_d_W;
    if (d_update_right_SV) {
        // The new d_W is the product of the current d_W extended by another row
        // and column and W.  The only new value in the extended version of d_W
        // that is non-zero is the new lower right value and it is 1.  We will
        // construct this product without explicitly forming the extended version of
        // d_W.
        new_d_W = new Matrix(d_num_rows_of_W+1, d_num_samples, false);
        for (int row = 0; row < d_num_rows_of_W; ++row) {
            for (int col = 0; col < d_num_samples; ++col) {
                double new_d_W_entry = 0.0;
                for (int entry = 0; entry < d_num_samples; ++entry) {
                    new_d_W_entry += d_W->item(row, entry)*W->item(entry, col);
                }
                new_d_W->item(row, col) = new_d_W_entry;
            }
        }
        for (int col = 0; col < d_num_samples; ++col) {
            new_d_W->item(d_num_rows_of_W, col) = W->item(d_num_samples, col);
        }
        delete d_W;
        d_W = new_d_W;
        ++d_num_rows_of_W;
    }

}

void
IncrementalSVDFastUpdate::addNewSample(
    const Vector* j,
    const Matrix* A,
    const Matrix* W,
    Matrix* sigma)
{
    CAROM_VERIFY(j != 0);
    CAROM_VERIFY(A != 0);
    CAROM_VERIFY(sigma != 0);

    // Add j as a new column of d_U.
    Matrix* newU = new Matrix(d_dim, d_num_samples+1, true);
    for (int row = 0; row < d_dim; ++row) {
        for (int col = 0; col < d_num_samples; ++col) {
            newU->item(row, col) = d_U->item(row, col);
        }
        newU->item(row, d_num_samples) = j->item(row);
    }
    delete d_U;
    d_U = newU;

    Matrix* new_d_W;
    if (d_update_right_SV) {
        // The new d_W is the product of the current d_W extended by another row
        // and column and W.  The only new value in the extended version of d_W
        // that is non-zero is the new lower right value and it is 1.  We will
        // construct this product without explicitly forming the extended version of
        // d_W.
        new_d_W = new Matrix(d_num_rows_of_W+1, d_num_samples+1, false);
        for (int row = 0; row < d_num_rows_of_W; ++row) {
            for (int col = 0; col < d_num_samples+1; ++col) {
                double new_d_W_entry = 0.0;
                for (int entry = 0; entry < d_num_samples; ++entry) {
                    new_d_W_entry += d_W->item(row, entry)*W->item(entry, col);
                }
                new_d_W->item(row, col) = new_d_W_entry;
            }
        }
        for (int col = 0; col < d_num_samples+1; ++col) {
            new_d_W->item(d_num_rows_of_W, col) = W->item(d_num_samples, col);
        }
        delete d_W;
        d_W = new_d_W;
    }

    // The new d_Up is the product of the current d_Up extended by another row
    // and column and A.  The only new value in the extended version of d_Up
    // that is non-zero is the new lower right value and it is 1.  We will
    // construct this product without explicitly forming the extended version of
    // d_Up.
    Matrix* new_d_Up = new Matrix(d_num_samples+1, d_num_samples+1, false);
    for (int row = 0; row < d_num_samples; ++row) {
        for (int col = 0; col < d_num_samples+1; ++col) {
            double new_d_Up_entry = 0.0;
            for (int entry = 0; entry < d_num_samples; ++entry) {
                new_d_Up_entry += d_Up->item(row, entry)*A->item(entry, col);
            }
            new_d_Up->item(row, col) = new_d_Up_entry;
        }
    }
    for (int col = 0; col < d_num_samples+1; ++col) {
        new_d_Up->item(d_num_samples, col) = A->item(d_num_samples, col);
    }
    delete d_Up;
    d_Up = new_d_Up;

    // d_S = sigma.
    delete d_S;
    int num_dim = std::min(sigma->numRows(), sigma->numColumns());
    d_S = new Vector(num_dim, false);
    for (int i = 0; i < num_dim; i++) {
        d_S->item(i) = sigma->item(i,i);
    }

    // We now have another sample.
    ++d_num_samples;
    ++d_num_rows_of_W;

    // Reorthogonalize if necessary.
    long int max_U_dim;
    if (d_num_samples > d_total_dim) {
        max_U_dim = d_num_samples;
    }
    else {
        max_U_dim = d_total_dim;
    }

}

}
