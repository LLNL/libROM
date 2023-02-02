/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Implementation of the MatrixInterpolator algorithm.

#include "MatrixInterpolator.h"

#include <limits.h>
#include <cmath>
#include "linalg/Matrix.h"
#include "linalg/scalapack_wrapper.h"
#include "mpi.h"

/* Use C++11 built-in shared pointers if available; else fallback to Boost. */
#if __cplusplus >= 201103L
#include <memory>
#else
#include <boost/shared_ptr.hpp>
#endif

/* Use automatically detected Fortran name-mangling scheme */
#define dposv CAROM_FC_GLOBAL(dposv, DPOSV)

extern "C" {
    // Solve a system of linear equations.
    void dposv(char*, int*, int*, double*, int*, double*, int*, int*);
}

using namespace std;

namespace CAROM {

MatrixInterpolator::MatrixInterpolator(std::vector<Vector*> parameter_points,
                                       std::vector<Matrix*> rotation_matrices,
                                       std::vector<Matrix*> reduced_matrices,
                                       int ref_point,
                                       std::string matrix_type,
                                       std::string rbf,
                                       std::string interp_method,
                                       double closest_rbf_val) :
    Interpolator(parameter_points,
                 rotation_matrices,
                 ref_point,
                 rbf,
                 interp_method,
                 closest_rbf_val)
{
    CAROM_VERIFY(reduced_matrices.size() == rotation_matrices.size());
    CAROM_VERIFY(matrix_type == "SPD" || matrix_type == "NS" || matrix_type == "R"
                 || matrix_type == "B");
    d_matrix_type = matrix_type;

    // Rotate the reduced matrices
    for (int i = 0; i < reduced_matrices.size(); i++)
    {
        // The ref_point does not need to be rotated
        if (i == d_ref_point)
        {
            d_rotated_reduced_matrices.push_back(reduced_matrices[i]);
            d_rotated_reduced_matrices_owned.push_back(false);
        }
        else
        {
            if (d_matrix_type == "B")
            {
                Matrix* AQ = reduced_matrices[i]->mult(rotation_matrices[i]);
                d_rotated_reduced_matrices.push_back(AQ);
            }
            else
            {
                Matrix* Q_tA = rotation_matrices[i]->transposeMult(reduced_matrices[i]);
                Matrix* Q_tAQ = Q_tA->mult(rotation_matrices[i]);
                delete Q_tA;
                d_rotated_reduced_matrices.push_back(Q_tAQ);
            }

            d_rotated_reduced_matrices_owned.push_back(true);
        }
    }
}

MatrixInterpolator::~MatrixInterpolator()
{
    CAROM_VERIFY(d_rotated_reduced_matrices.size() ==
                 d_rotated_reduced_matrices_owned.size());

    for (int i=0; i<d_rotated_reduced_matrices.size(); ++i)
    {
        if (d_rotated_reduced_matrices_owned[i])
            delete d_rotated_reduced_matrices[i];
    }

    for (auto matrix : d_gammas)
        delete matrix;

    delete d_x_half_power;
}

Matrix* MatrixInterpolator::interpolate(Vector* point, bool orthogonalize)
{
    Matrix* interpolated_matrix = NULL;
    if (d_matrix_type == "SPD")
    {
        interpolated_matrix = interpolateSPDMatrix(point);
    }
    else if (d_matrix_type == "NS")
    {
        interpolated_matrix = interpolateNonSingularMatrix(point);
    }
    else
    {
        interpolated_matrix = interpolateMatrix(point);
    }

    // Orthogonalize the interpolated matrix if requested.
    if (orthogonalize)
    {
        interpolated_matrix->orthogonalize();
    }
    return interpolated_matrix;
}

void MatrixInterpolator::obtainLambda()
{
    if (d_interp_method == "LS")
    {
        // Solving f = B*lambda
        Matrix* f_T = new Matrix(d_gammas[0]->numRows() * d_gammas[0]->numColumns(),
                                 d_gammas.size(), false);
        for (int i = 0; i < f_T->numRows(); i++)
        {
            for (int j = 0; j < f_T->numColumns(); j++)
            {
                f_T->item(i, j) = d_gammas[j]->getData()[i];
            }
        }

        // Obtain B matrix by calculating RBF.
        Matrix* B = new Matrix(d_gammas.size(), d_gammas.size(), false);
        for (int i = 0; i < B->numRows(); i++)
        {
            B->item(i, i) = 1.0;
            for (int j = i + 1; j < B->numColumns(); j++)
            {
                double res = obtainRBF(d_rbf, d_epsilon, d_parameter_points[i],
                                       d_parameter_points[j]);
                B->item(i, j) = res;
                B->item(j, i) = res;
            }
        }

        char uplo = 'U';
        int gamma_size = d_gammas.size();
        int num_elements = d_gammas[0]->numRows() * d_gammas[0]->numColumns();
        int info;

        dposv(&uplo, &gamma_size, &num_elements, B->getData(),  &gamma_size,
              f_T->getData(), &gamma_size, &info);
        if (info != 0)
        {
            std::cout << "Linear solve failed. Please choose a different epsilon value." <<
                      std::endl;
        }
        CAROM_VERIFY(info == 0);

        delete B;

        d_lambda_T = f_T;
    }
}

Matrix* MatrixInterpolator::obtainLogInterpolatedMatrix(
    std::vector<double>& rbf)
{
    Matrix* log_interpolated_matrix = new Matrix(
        d_rotated_reduced_matrices[d_ref_point]->numRows(),
        d_rotated_reduced_matrices[d_ref_point]->numColumns(),
        d_rotated_reduced_matrices[d_ref_point]->distributed());
    if (d_interp_method == "LS")
    {
        for (int i = 0; i < d_lambda_T->numRows(); i++)
        {
            for (int j = 0; j < rbf.size(); j++)
            {
                log_interpolated_matrix->getData()[i] += d_lambda_T->item(i, j) * rbf[j];
            }
        }
    }
    else if (d_interp_method == "IDW")
    {
        double sum = rbfWeightedSum(rbf);
        int num_elements = d_rotated_reduced_matrices[d_ref_point]->numRows() *
                           d_rotated_reduced_matrices[d_ref_point]->numColumns();
        for (int i = 0; i < num_elements; i++)
        {
            for (int j = 0; j < rbf.size(); j++)
            {
                log_interpolated_matrix->getData()[i] += d_gammas[j]->getData()[i] * rbf[j];
            }
            log_interpolated_matrix->getData()[i] /= sum;
        }
    }
    else if (d_interp_method == "LP")
    {
        int num_elements = d_rotated_reduced_matrices[d_ref_point]->numRows() *
                           d_rotated_reduced_matrices[d_ref_point]->numColumns();
        for (int i = 0; i < num_elements; i++)
        {
            for (int j = 0; j < rbf.size(); j++)
            {
                log_interpolated_matrix->getData()[i] += d_gammas[j]->getData()[i] * rbf[j];
            }
        }
    }
    return log_interpolated_matrix;
}

Matrix* MatrixInterpolator::interpolateSPDMatrix(Vector* point)
{
    if (d_gammas.size() == 0)
    {
        // Diagonalize X to work towards obtaining X^-1/2
        EigenPair ref_reduced_matrix_eigenpair = SymmetricRightEigenSolve(
                    d_rotated_reduced_matrices[d_ref_point]);
        Matrix* ref_reduced_matrix_ev = ref_reduced_matrix_eigenpair.ev;
        Matrix* ref_reduced_matrix_ev_inv = NULL;

        Matrix* ref_reduced_matrix_sqrt_eigs = new Matrix(
            ref_reduced_matrix_eigenpair.eigs.size(),
            ref_reduced_matrix_eigenpair.eigs.size(), false);
        for (int i = 0; i < ref_reduced_matrix_eigenpair.eigs.size(); i++)
        {
            ref_reduced_matrix_sqrt_eigs->item(i,
                                               i) = std::sqrt(ref_reduced_matrix_eigenpair.eigs[i]);
        }

        ref_reduced_matrix_ev->inverse(ref_reduced_matrix_ev_inv);

        // Obtain X^1/2
        Matrix* ref_reduced_matrix_ev_mult_sqrt_eig = ref_reduced_matrix_ev->mult(
                    ref_reduced_matrix_sqrt_eigs);
        d_x_half_power = ref_reduced_matrix_ev_mult_sqrt_eig->mult(
                             ref_reduced_matrix_ev_inv);
        Matrix* x_half_power_inv = NULL;

        // Obtain X^-1/2
        d_x_half_power->inverse(x_half_power_inv);

        delete ref_reduced_matrix_ev;
        delete ref_reduced_matrix_ev_inv;
        delete ref_reduced_matrix_sqrt_eigs;
        delete ref_reduced_matrix_ev_mult_sqrt_eig;

        // Obtain gammas for all points in the database.
        for (int i = 0; i < d_parameter_points.size(); i++)
        {
            // For the ref point, gamma is the zero matrix
            if (i == d_ref_point)
            {
                Matrix* gamma = new Matrix(x_half_power_inv->numRows(),
                                           x_half_power_inv->numColumns(), x_half_power_inv->distributed());
                d_gammas.push_back(gamma);
            }
            else
            {
                Matrix* x_half_power_inv_mult_y = x_half_power_inv->mult(
                                                      d_rotated_reduced_matrices[i]);

                // Obtain X^-1/2*Y*X^-1/2
                Matrix* x_half_power_inv_mult_y_mult_x_half_power_inv =
                    x_half_power_inv_mult_y->mult(x_half_power_inv);
                delete x_half_power_inv_mult_y;

                // Diagonalize X^-1/2*Y*X^-1/2 to obtain the log of this matrix.
                // Diagonalize YX^-1 to obtain log of this matrix.
                // Following https://en.wikipedia.org/wiki/Logarithm_of_a_matrix
                // 1. Diagonalize M to obtain M' = V^-1*M*V. M' are the eigenvalues
                // of M and V are the eigenvectors of M.
                // 2. log M = V(log M')V^-1
                EigenPair log_eigenpair = SymmetricRightEigenSolve(
                                              x_half_power_inv_mult_y_mult_x_half_power_inv);
                delete x_half_power_inv_mult_y_mult_x_half_power_inv;
                Matrix* log_ev = log_eigenpair.ev;
                Matrix* log_ev_inv = NULL;

                Matrix* log_eigs = new Matrix(log_eigenpair.eigs.size(),
                                              log_eigenpair.eigs.size(), false);
                for (int i = 0; i < log_eigenpair.eigs.size(); i++)
                {
                    if (log_eigenpair.eigs[i] < 0)
                    {
                        if (d_rank == 0) std::cout << "Some eigenvalues of this matrix are negative, "
                                                       <<
                                                       "which leads to NaN values when taking the log. Aborting." << std::endl;
                        CAROM_VERIFY(log_eigenpair.eigs[i] > 0);
                    }
                    log_eigs->item(i, i) = std::log(log_eigenpair.eigs[i]);
                }

                // Invert matrix.
                log_ev->inverse(log_ev_inv);

                // Perform log mapping.
                Matrix* log_ev_mult_log_eig = log_ev->mult(log_eigs);
                Matrix* gamma = log_ev_mult_log_eig->mult(log_ev_inv);
                d_gammas.push_back(gamma);

                delete log_ev;
                delete log_ev_inv;
                delete log_eigs;
                delete log_ev_mult_log_eig;
            }
        }

        delete x_half_power_inv;

        // Obtain lambda for the P interpolation matrix
        obtainLambda();
    }

    // Obtain distances from database points to new point
    std::vector<double> rbf = obtainRBFToTrainingPoints(d_parameter_points,
                              d_interp_method,
                              d_rbf, d_epsilon, point);

    // Interpolate gammas to get gamma for new point
    Matrix* log_interpolated_matrix = obtainLogInterpolatedMatrix(rbf);

    // Diagonalize the new gamma so we can exponentiate it
    // Diagonalize X to obtain exp(X) of this matrix.
    // Following https://en.wikipedia.org/wiki/Matrix_exponential
    // 1. Diagonalize M to obtain M' = V^-1*M*V. M' are the eigenvalues
    // of M and V are the eigenvectors of M.
    // 2. exp M = V(exp M')V^-1
    EigenPair exp_eigenpair = SymmetricRightEigenSolve(log_interpolated_matrix);
    delete log_interpolated_matrix;

    Matrix* exp_ev = exp_eigenpair.ev;
    Matrix* exp_ev_inv = NULL;

    Matrix* exp_eigs = new Matrix(exp_eigenpair.eigs.size(),
                                  exp_eigenpair.eigs.size(), false);
    for (int i = 0; i < exp_eigenpair.eigs.size(); i++)
    {
        exp_eigs->item(i, i) = std::exp(exp_eigenpair.eigs[i]);
    }

    exp_ev->inverse(exp_ev_inv);
    Matrix* exp_ev_mult_exp_eig = exp_ev->mult(exp_eigs);

    // Exponentiate gamma
    Matrix* exp_gamma = exp_ev_mult_exp_eig->mult(exp_ev_inv);

    delete exp_ev;
    delete exp_ev_inv;
    delete exp_eigs;
    delete exp_ev_mult_exp_eig;

    // Obtain exp mapping by doing X^1/2*exp(gamma)*X^1/2
    Matrix* x_half_power_mult_exp_gamma = d_x_half_power->mult(exp_gamma);
    Matrix* interpolated_matrix = x_half_power_mult_exp_gamma->mult(d_x_half_power);

    delete  x_half_power_mult_exp_gamma;
    return interpolated_matrix;

}

Matrix* MatrixInterpolator::interpolateNonSingularMatrix(Vector* point)
{
    if (d_gammas.size() == 0)
    {
        // Invert X
        Matrix* ref_matrix_inv = NULL;
        d_rotated_reduced_matrices[d_ref_point]->inverse(ref_matrix_inv);

        for (int i = 0; i < d_parameter_points.size(); i++)
        {
            // For ref_point, gamma is the zero matrix
            if (i == d_ref_point)
            {
                Matrix* gamma = new Matrix(ref_matrix_inv->numRows(),
                                           ref_matrix_inv->numColumns(), ref_matrix_inv->distributed());
                d_gammas.push_back(gamma);
            }
            else
            {
                Matrix* y_mult_ref_matrix_inv = d_rotated_reduced_matrices[i]->mult(
                                                    ref_matrix_inv);

                // Diagonalize YX^-1 to obtain log of this matrix.
                // Following https://en.wikipedia.org/wiki/Logarithm_of_a_matrix
                // 1. Diagonalize M to obtain M' = V^-1*M*V. M' are the eigenvalues
                // of M and V are the eigenvectors of M.
                // 2. log M = V(log M')V^-1
                EigenPair log_eigenpair = SymmetricRightEigenSolve(y_mult_ref_matrix_inv);
                delete y_mult_ref_matrix_inv;
                Matrix* log_ev = log_eigenpair.ev;
                Matrix* log_ev_inv = NULL;

                Matrix* log_eigs = new Matrix(log_eigenpair.eigs.size(),
                                              log_eigenpair.eigs.size(), false);
                for (int i = 0; i < log_eigenpair.eigs.size(); i++)
                {
                    if (log_eigenpair.eigs[i] < 0)
                    {
                        if (d_rank == 0) std::cout << "Some eigenvalues of this matrix are " <<
                                                       "negative, which leads to NaN values when taking the log. Aborting." <<
                                                       std::endl;
                        CAROM_VERIFY(log_eigenpair.eigs[i] > 0);
                    }
                    log_eigs->item(i, i) = std::log(log_eigenpair.eigs[i]);
                }

                log_ev->inverse(log_ev_inv);

                // Perform log mapping.
                Matrix* log_ev_mult_log_eig = log_ev->mult(log_eigs);
                Matrix* gamma = log_ev_mult_log_eig->mult(log_ev_inv);
                d_gammas.push_back(gamma);

                delete log_ev;
                delete log_ev_inv;
                delete log_eigs;
                delete log_ev_mult_log_eig;
            }
        }

        delete ref_matrix_inv;

        // Obtain lambda for the P interpolation matrix
        obtainLambda();
    }

    // Obtain distances from database points to new point
    std::vector<double> rbf = obtainRBFToTrainingPoints(d_parameter_points,
                              d_interp_method, d_rbf, d_epsilon, point);

    // Interpolate gammas to get gamma for new point
    Matrix* log_interpolated_matrix = obtainLogInterpolatedMatrix(rbf);

    // Diagonalize the new gamma so we can exponentiate it
    // Diagonalize X to obtain exp(X) of this matrix.
    // Following https://en.wikipedia.org/wiki/Matrix_exponential
    // 1. Diagonalize M to obtain M' = V^-1*M*V. M' are the eigenvalues
    // of M and V are the eigenvectors of M.
    // 2. exp M = V(exp M')V^-1
    EigenPair exp_eigenpair = SymmetricRightEigenSolve(log_interpolated_matrix);
    delete log_interpolated_matrix;
    Matrix* exp_ev = exp_eigenpair.ev;
    Matrix* exp_ev_inv = NULL;

    Matrix* exp_eigs = new Matrix(exp_eigenpair.eigs.size(),
                                  exp_eigenpair.eigs.size(), false);
    for (int i = 0; i < exp_eigenpair.eigs.size(); i++)
    {
        exp_eigs->item(i, i) = std::exp(exp_eigenpair.eigs[i]);
    }

    // Invert matrix.
    exp_ev->inverse(exp_ev_inv);

    // Perform log mapping.
    Matrix* exp_ev_mult_exp_eig = exp_ev->mult(exp_eigs);

    // Exponentiate gamma
    Matrix* exp_gamma = exp_ev_mult_exp_eig->mult(exp_ev_inv);

    delete exp_ev;
    delete exp_ev_inv;
    delete exp_eigs;
    delete exp_ev_mult_exp_eig;

    // Obtain exp mapping by doing exp(gamma)*X
    Matrix* interpolated_matrix = exp_gamma->mult(
                                      d_rotated_reduced_matrices[d_ref_point]);
    delete exp_gamma;

    return interpolated_matrix;
}

Matrix* MatrixInterpolator::interpolateMatrix(Vector* point)
{
    if (d_gammas.size() == 0)
    {
        for (int i = 0; i < d_parameter_points.size(); i++)
        {
            // For ref point, gamma is the zero matrix.
            if (i == d_ref_point)
            {
                Matrix* gamma = new Matrix(d_rotated_reduced_matrices[d_ref_point]->numRows(),
                                           d_rotated_reduced_matrices[d_ref_point]->numColumns(),
                                           d_rotated_reduced_matrices[d_ref_point]->distributed());
                d_gammas.push_back(gamma);
            }
            else
            {
                // Gamma is Y - X
                Matrix* gamma = new Matrix(*d_rotated_reduced_matrices[i]);
                *gamma -= *d_rotated_reduced_matrices[d_ref_point];
                d_gammas.push_back(gamma);
            }
        }

        // Obtain lambda for the P interpolation matrix
        obtainLambda();
    }
    // Obtain distances from database points to new point
    std::vector<double> rbf = obtainRBFToTrainingPoints(d_parameter_points,
                              d_interp_method,
                              d_rbf, d_epsilon, point);

    // Interpolate gammas to get gamma for new point
    Matrix* interpolated_matrix = obtainLogInterpolatedMatrix(rbf);

    // The exp mapping is X + the interpolated gamma
    *interpolated_matrix += *d_rotated_reduced_matrices[d_ref_point];
    return interpolated_matrix;
}

}
