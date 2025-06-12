/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Implementation of the Interpolator algorithm.

#include "Interpolator.h"

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

using namespace std;

namespace CAROM {

Interpolator::Interpolator(const std::vector<Vector> & parameter_points,
                           const std::vector<std::shared_ptr<Matrix>> & rotation_matrices,
                           int ref_point,
                           std::string rbf,
                           std::string interp_method,
                           double closest_rbf_val,
                           bool compute_gradients)
{
    CAROM_VERIFY(parameter_points.size() == rotation_matrices.size());
    CAROM_VERIFY(parameter_points.size() > 1);
    CAROM_VERIFY(rbf == "G" || rbf == "IQ" || rbf == "IMQ");
    CAROM_VERIFY(interp_method == "LS" || interp_method == "IDW"
                 || interp_method == "LP");
    CAROM_VERIFY(closest_rbf_val >= 0.0 && closest_rbf_val <= 1.0);

    // Get the rank of this process, and the number of processors.
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

    d_compute_gradients = compute_gradients;
    d_parameter_points = parameter_points;
    d_rotation_matrices = rotation_matrices;
    d_ref_point = ref_point;
    d_rbf = rbf;
    d_interp_method = interp_method;
    d_epsilon = convertClosestRBFToEpsilon(parameter_points, rbf, closest_rbf_val);
}

std::vector<double> obtainRBFToTrainingPoints(
    const std::vector<Vector> & parameter_points,
    const std::string & interp_method, const std::string & rbf, double epsilon,
    const Vector & point)
{
    std::vector<double> rbfs;
    if (interp_method == "LS")
    {
        for (int i = 0; i < parameter_points.size(); i++)
        {
            rbfs.push_back(obtainRBF(rbf, epsilon, point, parameter_points[i]));
        }
    }
    else if (interp_method == "IDW")
    {
        bool distance_is_zero = false;
        for (int i = 0; i < parameter_points.size(); i++)
        {
            rbfs.push_back(obtainRBF(rbf, epsilon, point, parameter_points[i]));
            if (rbfs.back() == 1.0)
            {
                distance_is_zero = true;
            }
        }
        if (distance_is_zero)
        {
            for (int i = 0; i < rbfs.size(); i++)
            {
                if (rbfs[i] != 1.0)
                {
                    rbfs[i] = 0.0;
                }
            }
        }
    }
    else if (interp_method == "LP")
    {
        for (int i = 0; i < parameter_points.size(); i++)
        {
            double coeff;
            bool first = true;
            for (int j = 0; j < parameter_points.size(); j++)
            {
                if (i == j)
                {
                    continue;
                }

                Vector numerator_vec, denomenator_vec;
                point.minus(parameter_points[j], numerator_vec);
                parameter_points[i].minus(parameter_points[j], denomenator_vec);
                double numerator = numerator_vec.norm();
                double denomenator = denomenator_vec.norm();

                if (first)
                {
                    coeff = (numerator / denomenator);
                    first = false;
                }
                else
                {
                    coeff *= (numerator / denomenator);
                }
            }
            rbfs.push_back(coeff);
        }
    }
    return rbfs;
}

std::vector<double> obtainRBFGradientToTrainingPoints(
    const std::vector<Vector> & parameter_points,
    const std::string & interp_method, const std::string & rbf, double epsilon,
    const Vector & point, const int index)
{
    std::vector<double> rbfs;
    if (interp_method == "LS")
    {
        for (int i = 0; i < parameter_points.size(); i++)
        {
            rbfs.push_back(obtainRBFGradient(rbf, epsilon, point, parameter_points[i],
                                             index));
        }
    }
    else
    {
        CAROM_ERROR("Interpolated gradients are only implemented for \"LS\"");
    }
    return rbfs;
}

double obtainRBFGradient(std::string rbf, double epsilon, const Vector & point1,
                         const Vector & point2, const int index)
{
    Vector diff;
    point1.minus(point2, diff);
    const double eps_norm_squared = epsilon * epsilon * diff.norm2();
    double res = 0.0;

    // Gaussian RBF
    if (rbf == "G")
    {
        //Derivative of Gaussian RBF, res = std::exp(-eps_norm_squared);
        res = -2.0*epsilon*epsilon*diff.item(index)*std::exp(-eps_norm_squared);
    }
    // Inverse quadratic RBF
    else if (rbf == "IQ")
    {
        //Derivative of Inverse Quadratic RBF, res = 1.0 / (1.0 + eps_norm_squared);
        res = -2.0*epsilon*epsilon*diff.item(index)/((1.0+eps_norm_squared)*
                (1.0+eps_norm_squared));

    }
    // Inverse multiquadric RBF
    else if (rbf == "IMQ")
    {
        //Derivative of Inverse multiquadritic RBF, res = 1.0 / std::sqrt(1.0 + eps_norm_squared);
        res = -epsilon*epsilon*diff.item(index)/(std::sqrt(1.0+eps_norm_squared)*
                (1.0+eps_norm_squared));
    }

    return res;
}

double rbfWeightedSum(std::vector<double>& rbf)
{
    double sum = 0.0;
    for (int i = 0; i < rbf.size(); i++)
    {
        sum += rbf[i];
    }
    return sum;
}

double obtainRBF(std::string rbf, double epsilon, const Vector & point1,
                 const Vector & point2)
{
    Vector diff;
    point1.minus(point2, diff);
    const double eps_norm_squared = epsilon * epsilon * diff.norm2();
    double res = 0.0;

    // Gaussian RBF
    if (rbf == "G")
    {
        res = std::exp(-eps_norm_squared);
    }
    // Inverse quadratic RBF
    else if (rbf == "IQ")
    {
        res = 1.0 / (1.0 + eps_norm_squared);
    }
    // Inverse multiquadric RBF
    else if (rbf == "IMQ")
    {
        res = 1.0 / std::sqrt(1.0 + eps_norm_squared);
    }

    return res;
}

double convertClosestRBFToEpsilon(const std::vector<Vector> & parameter_points,
                                  const std::string & rbf, double closest_rbf_val)
{
    double closest_point_dist = INT_MAX;
    double epsilon;
    for (int i = 0; i < parameter_points.size(); i++)
    {
        for (int j = 0; j < parameter_points.size(); j++)
        {
            if (i == j)
            {
                continue;
            }

            Vector diff;
            parameter_points[i].minus(parameter_points[j], diff);
            double dist = diff.norm2();
            if (dist < closest_point_dist)
            {
                closest_point_dist = dist;

                // Gaussian RBF
                if (rbf == "G")
                {
                    epsilon = std::sqrt(-std::log(closest_rbf_val) / diff.norm2());
                }
                // Inverse quadratic RBF
                else if (rbf == "IQ")
                {
                    epsilon = std::sqrt(((1.0 / closest_rbf_val) - 1.0) / diff.norm2());
                }
                // Inverse multiquadric RBF
                else if (rbf == "IMQ")
                {
                    epsilon = std::sqrt((std::pow(1.0 / closest_rbf_val, 2) - 1.0) / diff.norm2());
                }
            }
        }
    }

    return epsilon;
}

std::vector<std::shared_ptr<Matrix>> obtainRotationMatrices(
                                      const std::vector<Vector> & parameter_points,
                                      const std::vector<std::shared_ptr<Matrix>> & bases,
                                      int ref_point)
{
    // Get the rank of this process, and the number of processors.
    int mpi_init;
    int rank;
    int num_procs;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    std::vector<std::shared_ptr<Matrix>> rotation_matrices;

    // Obtain the rotation matrices to rotate the bases into
    // the same generalized coordinate space.
    for (int i = 0; i < parameter_points.size(); i++)
    {
        CAROM_VERIFY(bases[i]->numRows() == bases[ref_point]->numRows());
        CAROM_VERIFY(bases[i]->numColumns() == bases[ref_point]->numColumns());
        CAROM_VERIFY(bases[i]->distributed() == bases[ref_point]->distributed());

        // If at ref point, the rotation_matrix is the identity matrix
        // since the ref point doesn't need to be rotated.
        if (i == ref_point)
        {
            std::shared_ptr<Matrix> identity_matrix(new Matrix(bases[i]->numColumns(),
                                                    bases[i]->numColumns(), false));
            for (int j = 0; j < identity_matrix->numColumns(); j++) {
                identity_matrix->item(j, j) = 1.0;
            }
            rotation_matrices.push_back(identity_matrix);
            continue;
        }

        std::unique_ptr<Matrix> basis_mult_basis = bases[i]->transposeMult(
                    *bases[ref_point]);
        Matrix basis(basis_mult_basis->numRows(),
                     basis_mult_basis->numColumns(), false);
        Matrix basis_right(basis_mult_basis->numColumns(),
                           basis_mult_basis->numColumns(), false);

        // We need to compute the SVD of basis_mult_basis. Since it is
        // undistributed due to the transposeMult, let's use lapack's serial SVD
        // on rank 0 only.
        if (rank == 0)
        {
            Vector sv(basis_mult_basis->numColumns(), false);
            SerialSVD(basis_mult_basis.get(), &basis_right, &sv, &basis);
        }

        // Broadcast U and V which are computed only on root.
        MPI_Bcast(basis.getData(), basis.numRows() * basis.numColumns(), MPI_DOUBLE,
                  0, MPI_COMM_WORLD);
        MPI_Bcast(basis_right.getData(),
                  basis_right.numRows() * basis_right.numColumns(), MPI_DOUBLE, 0,
                  MPI_COMM_WORLD);

        // Obtain the rotation matrix.
        rotation_matrices.push_back(basis.mult(basis_right));
    }

    return rotation_matrices;
}

std::unique_ptr<std::vector<Vector>>
                                  scalarsToVectors(const std::vector<double> & s)
{
    std::vector<Vector> *v = new std::vector<Vector>(s.size());
    for (int i=0; i<s.size(); ++i)
    {
        (*v)[i].setSize(1);
        (*v)[i](0) = s[i];
    }

    return std::unique_ptr<std::vector<Vector>>(v);
}

}
