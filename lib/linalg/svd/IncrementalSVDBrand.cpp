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

#include "IncrementalSVDBrand.h"
#include "utils/HDFDatabase.h"

#include "mpi.h"

#include <cmath>
#include <limits>

namespace CAROM {

IncrementalSVDBrand::IncrementalSVDBrand(
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

IncrementalSVDBrand::~IncrementalSVDBrand()
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

const Matrix*
IncrementalSVDBrand::getSpatialBasis()
{
    updateSpatialBasis(); // WARNING: this is costly

    CAROM_ASSERT(d_basis != 0);
    return d_basis;
}

const Matrix*
IncrementalSVDBrand::getTemporalBasis()
{
    updateTemporalBasis();

    CAROM_ASSERT(d_basis_right != 0);
    return d_basis_right;
}

void
IncrementalSVDBrand::buildInitialSVD(
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

}

bool
IncrementalSVDBrand::buildIncrementalSVD(
    double* u, bool add_without_increase)
{
    CAROM_VERIFY(u != 0);

    // Compute the projection error
    // (accurate down to the machine precision)
    Vector u_vec(u, d_dim, true);
    Vector e_proj(u, d_dim, true);
    e_proj -= *(d_U->mult(d_U->transposeMult(e_proj))); // Gram-Schmidt
    e_proj -= *(d_U->mult(d_U->transposeMult(e_proj))); // Re-orthogonalization

    double k = e_proj.inner_product(e_proj);
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
    Vector* U_mult_u = new Vector(d_U->transposeMult(u_vec)->getData(),
                                  d_num_samples,
                                  false);
    Vector* l = d_Up->transposeMult(U_mult_u);
    constructQ(Q, l, k);
    delete U_mult_u, l;

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
            Vector* j = new Vector(e_proj.getData(), d_dim, false);
            for (int i = 0; i < d_dim; ++i) {
                j->item(i) /= k;
            }

            // addNewSample will assign sigma to d_S hence it should not be
            // deleted upon return.
            addNewSample(j, A, W, sigma);
            delete j;
        }
        delete A;
        delete W;
    }
    else {
        delete A;
        delete W;
        delete sigma;
    }
    return result;
}

void
IncrementalSVDBrand::updateSpatialBasis()
{
    d_basis = d_U->mult(d_Up);

    // remove the smallest singular value if it is smaller than d_singular_value_tol
    if ( (d_singular_value_tol != 0.0) &&
            (d_S->item(d_num_samples-1) < d_singular_value_tol) &&
            (d_num_samples != 1) ) {

        if (d_rank == 0) {
            std::cout <<
                      "removing a spatial basis corresponding to the small singular value!\n";
        }

        Matrix* d_basis_new = new Matrix(d_dim, d_num_samples-1,
                                         d_basis->distributed());
        for (int row = 0; row < d_dim; ++row) {
            for (int col = 0; col < d_num_samples-1; ++col) {
                d_basis_new->item(row, col) = d_basis->item(row,col);
            }
        }
        delete d_basis;
        d_basis = d_basis_new;
    }

    // Reorthogonalize if necessary.
    // (not likely to be called anymore but left for safety)
    if (fabs(checkOrthogonality(d_basis)) >
            std::numeric_limits<double>::epsilon()*static_cast<double>(d_num_samples)) {
        d_basis->orthogonalize();
    }

}

void
IncrementalSVDBrand::updateTemporalBasis()
{
    delete d_basis_right;
    d_basis_right = new Matrix(*d_W);

    // remove the smallest singular value if it is smaller than d_singular_value_tol
    if ( (d_singular_value_tol != 0.0) &&
            (d_S->item(d_num_samples-1) < d_singular_value_tol) &&
            (d_num_samples != 1) ) {

        if (d_rank == 0) {
            std::cout <<
                      "removing a temporal basis corresponding to the small singular value!\n";
        }

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

    // Reorthogonalize if necessary.
    // (not likely to be called anymore but left for safety)
    if (fabs(checkOrthogonality(d_basis_right)) >
            std::numeric_limits<double>::epsilon()*d_num_samples) {
        d_basis_right->orthogonalize();
    }

}

void
IncrementalSVDBrand::computeBasis()
{
    if(d_rank == 0) {
        std::cout << "d_num_samples = " << d_num_samples << "\n";
        std::cout << "d_num_rows_of_W = " << d_num_rows_of_W << "\n";
        std::cout << "d_singular_value_tol = " << d_singular_value_tol << "\n";
        std::cout << "smallest SV = " << d_S->item(d_num_samples-1) << "\n";
        if (d_num_samples > 1) {
            std::cout << "next smallest SV = " << d_S->item(d_num_samples-2) << "\n";
        }
    }

    updateSpatialBasis();
    if (d_update_right_SV)
    {
        updateTemporalBasis();
    }

    // remove the smallest singular value if it is smaller than d_singular_value_tol
    if ( (d_singular_value_tol != 0.0) &&
            (d_S->item(d_num_samples-1) < d_singular_value_tol) &&
            (d_num_samples != 1) ) {

        --d_num_samples;
    }
}

void
IncrementalSVDBrand::addLinearlyDependentSample(
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
IncrementalSVDBrand::addNewSample(
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
    if (fabs(checkOrthogonality(d_Up)) >
            std::numeric_limits<double>::epsilon()*static_cast<double>(max_U_dim)) {
        d_Up->orthogonalize();
    }
    if (fabs(checkOrthogonality(d_U)) >
            std::numeric_limits<double>::epsilon()*static_cast<double>(max_U_dim)) {
        d_U->orthogonalize(); // Will not be called, but just in case
    }

    if(d_update_right_SV)
    {
        if (fabs(checkOrthogonality(d_W)) >
                std::numeric_limits<double>::epsilon()*d_num_samples) {
            d_W->orthogonalize();
        }
    }

}

}
