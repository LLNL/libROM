/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Implementation of the VectorInterpolator algorithm.

#include "VectorInterpolator.h"

#include <limits.h>
#include <cmath>
#include "linalg/Matrix.h"
#include "linalg/Vector.h"
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

VectorInterpolator::VectorInterpolator(std::vector<Vector*> parameter_points,
                                       std::vector<Matrix*> rotation_matrices,
                                       std::vector<Vector*> reduced_vectors,
                                       int ref_point,
                                       std::string rbf,
                                       std::string interp_method,
                                       double epsilon) :
    Interpolator(parameter_points,
                 rotation_matrices,
                 ref_point,
                 rbf,
                 interp_method,
                 epsilon)
{
    CAROM_VERIFY(reduced_vectors.size() == rotation_matrices.size());

    // Rotate the reduced vectors
    for (int i = 0; i < reduced_vectors.size(); i++)
    {
        // The ref_point does not need to be rotated
        if (i == d_ref_point)
        {
            d_rotated_reduced_vectors.push_back(reduced_vectors[i]);
        }
        else
        {
            Vector* Q_tA = rotation_matrices[i]->transposeMult(reduced_vectors[i]);
            d_rotated_reduced_vectors.push_back(Q_tA);
        }
    }
}

void VectorInterpolator::obtainLambda()
{
    if (d_interp_method == "LS")
    {
        // Solving f = B*lambda
        Matrix* f_T = new Matrix(d_gammas[0]->dim(), d_gammas.size(), false);
        for (int i = 0; i < f_T->numRows(); i++)
        {
            for (int j = 0; j < f_T->numColumns(); j++)
            {
                f_T->item(i, j) = d_gammas[j]->getData()[i];
            }
        }

        // Obtain B vector by calculating RBF.
        Matrix* B = new Matrix(d_gammas.size(), d_gammas.size(), false);
        for (int i = 0; i < B->numRows(); i++)
        {
            B->item(i, i) = 1.0;
            for (int j = i + 1; j < B->numColumns(); j++)
            {
                double res = obtainRBF(d_parameter_points[i], d_parameter_points[j]);
                B->item(i, j) = res;
                B->item(j, i) = res;
            }
        }

        char uplo = 'U';
        int gamma_size = d_gammas.size();
        int num_elements = d_gammas[0]->dim();
        int info;

        dposv(&uplo, &gamma_size, &num_elements, B->getData(),  &gamma_size, f_T->getData(), &gamma_size, &info);
        CAROM_VERIFY(info == 0);

        delete B;

        d_lambda_T = f_T;
    }
}

Vector* VectorInterpolator::obtainLogInterpolatedVector(std::vector<double> rbf)
{
    Vector* log_interpolated_vector = new Vector(d_rotated_reduced_vectors[d_ref_point]->dim(), d_rotated_reduced_vectors[d_ref_point]->distributed());
    if (d_interp_method == "LS")
    {
        for (int i = 0; i < d_lambda_T->numRows(); i++)
        {
            for (int j = 0; j < rbf.size(); j++)
            {
                log_interpolated_vector->getData()[i] += d_lambda_T->item(i, j) * rbf[j];
            }
        }
    }
    else if (d_interp_method == "IDW")
    {
        double sum = rbfWeightedSum(rbf);
        for (int i = 0; i < d_rotated_reduced_vectors[d_ref_point]->dim(); i++)
        {
            for (int j = 0; j < rbf.size(); j++)
            {
                log_interpolated_vector->getData()[i] += d_gammas[j]->getData()[i] * rbf[j];
            }
            log_interpolated_vector->getData()[i] /= sum;
        }
    }
    else if (d_interp_method == "LP")
    {
        for (int i = 0; i < d_rotated_reduced_vectors[d_ref_point]->dim(); i++)
        {
            for (int j = 0; j < rbf.size(); j++)
            {
                log_interpolated_vector->getData()[i] += d_gammas[j]->getData()[i] * rbf[j];
            }
        }
    }

    return log_interpolated_vector;
}

Vector* VectorInterpolator::interpolate(Vector* point)
{
    if (d_gammas.size() == 0)
    {

        for (int i = 0; i < d_parameter_points.size(); i++)
        {
            // For ref point, gamma is the zero vector.
            if (i == d_ref_point)
            {
                Vector* gamma = new Vector(d_rotated_reduced_vectors[d_ref_point]->dim(), d_rotated_reduced_vectors[d_ref_point]->distributed());
                d_gammas.push_back(gamma);
            }
            else
            {
                // Gamma is Y - X
                Vector* gamma = d_rotated_reduced_vectors[i]->minus(*d_rotated_reduced_vectors[d_ref_point]);
                d_gammas.push_back(gamma);
            }
        }

        // Obtain lambda for the P interpolation vector
        obtainLambda();
    }

    // Obtain distances from database points to new point
    std::vector<double> rbf = obtainRBFToTrainingPoints(point);

    // Interpolate gammas to get gamma for new point
    Vector* log_interpolated_vector = obtainLogInterpolatedVector(rbf);

    // The exp mapping is X + the interpolated gamma
    Vector* interpolated_vector = d_rotated_reduced_vectors[d_ref_point]->plus(log_interpolated_vector);
    delete log_interpolated_vector;
    return interpolated_vector;
}

}
