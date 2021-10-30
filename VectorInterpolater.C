/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Implementation of the VectorInterpolater algorithm.

#include "VectorInterpolater.h"

#include <limits.h>
#include <cmath>
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

/* Use automatically detected Fortran name-mangling scheme */
#define dgesv CAROM_FC_GLOBAL(dgesv, DGESV)

extern "C" {
    // Solve a system of linear equations.
    void dgesv(int*, int*, double*, int*, int*, double*, int*, int*);
}

using namespace std;

namespace CAROM {

VectorInterpolater::VectorInterpolater(std::vector<Vector*> parameter_points,
                                       std::vector<Matrix*> rotation_matrices,
                                       std::vector<Vector*> reduced_vectors,
                                       int ref_point,
                                       double epsilon) :
    Interpolater(parameter_points,
                 rotation_matrices,
                 ref_point,
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

void VectorInterpolater::obtainLambda(std::vector<Vector*> gammas)
{

    // Solving f = B*lambda
    Matrix* f_T = new Matrix(gammas[0]->dim(), gammas.size(), false);
    for (int i = 0; i < f_T->numRows(); i++)
    {
        for (int j = 0; j < f_T->numColumns(); j++)
        {
            f_T->item(i, j) = gammas[j]->getData()[i];
        }
    }

    // Obtain B vector by calculating inverse quadratic distance.
    Matrix* B = new Matrix(gammas.size(), gammas.size(), false);
    for (int i = 0; i < B->numRows(); i++)
    {
        for (int j = i + 1; j < B->numColumns(); j++)
        {
            Vector diff;
            d_parameter_points[i]->minus(*d_parameter_points[j], diff);
            double res = 1.0 / (1.0 + std::pow(diff.norm(), 2));
            B->item(i, j) = res;
            B->item(j, i) = res;
        }
    }

    int gamma_size = gammas.size();
    int num_elements = gammas[0]->dim();
    int* ipiv = new int[gamma_size];
    int info;

    dgesv(&gamma_size, &num_elements, B->getData(),  &gamma_size, ipiv, f_T->getData(), &gamma_size, &info);
    CAROM_VERIFY(info == 0);

    delete [] ipiv;
    delete B;

    d_lambda_T = f_T;
}

Vector* VectorInterpolater::obtainLogInterpolatedVector(std::vector<double> inv_q)
{
    Vector* log_interpolated_vector = new Vector(d_rotated_reduced_vectors[d_ref_point]->dim(), d_rotated_reduced_vectors[d_ref_point]->distributed());
    CAROM_VERIFY(d_rotated_reduced_vectors[d_ref_point]->dim() == d_lambda_T->numRows());
    for (int i = 0; i < d_lambda_T->numRows(); i++)
    {
        for (int j = 0; j < inv_q.size(); j++)
        {
            log_interpolated_vector->getData()[i] += d_lambda_T->item(i, j) * inv_q[j];
        }
    }

    return log_interpolated_vector;
}

Vector* VectorInterpolater::interpolate(Vector* point)
{
    if (d_lambda_T == NULL)
    {
        // Perform log mapping by doing Y - X.
        std::vector<Vector*> gammas;

        for (int i = 0; i < d_parameter_points.size(); i++)
        {
            // For ref point, gamma is the zero vector.
            if (i == d_ref_point)
            {
                Vector* gamma = new Vector(d_rotated_reduced_vectors[d_ref_point]->dim(), d_rotated_reduced_vectors[d_ref_point]->distributed());
                gammas.push_back(gamma);
            }
            else
            {
                // Gamma is Y - X
                Vector* gamma = d_rotated_reduced_vectors[i]->minus(*d_rotated_reduced_vectors[d_ref_point]);
                gammas.push_back(gamma);
            }
        }

        // Obtain lambda for the P interpolation vector
        obtainLambda(gammas);
    }

    // Obtain distances from database points to new point
    std::vector<double> inv_q = obtainRBF(point);

    // Interpolate gammas to get gamma for new point
    Vector* log_interpolated_vector = obtainLogInterpolatedVector(inv_q);

    // The exp mapping is X + the interpolated gamma
    Vector* interpolated_vector = d_rotated_reduced_vectors[d_ref_point]->plus(log_interpolated_vector);
    delete log_interpolated_vector;
    return interpolated_vector;
}

}
