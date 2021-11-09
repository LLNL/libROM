/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
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
#include "Matrix.h"
#include "scalapack_wrapper.h"
#include "mpi.h"

/* Use C++11 built-in shared pointers if available; else fallback to Boost. */
#if __cplusplus >= 201103L
#include <memory>
#else
#include <boost/shared_ptr.hpp>
#endif

using namespace std;

namespace CAROM {

Interpolator::Interpolator(std::vector<Vector*> parameter_points,
                           std::vector<Matrix*> rotation_matrices,
                           int ref_point,
                           std::string rbf,
                           double epsilon)
{
    CAROM_VERIFY(parameter_points.size() == rotation_matrices.size());
    CAROM_VERIFY(parameter_points.size() > 0);
    CAROM_VERIFY(rbf == "G" || rbf == "IQ" || rbf == "MQ" || rbf == "IMQ");

    // Get the rank of this process, and the number of processors.
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

    d_parameter_points = parameter_points;
    d_rotation_matrices = rotation_matrices;
    d_ref_point = ref_point;
    d_lambda_T = NULL;
    d_rbf = rbf;
    d_epsilon = epsilon;
}

std::vector<double> Interpolator::obtainRBFToTrainingPoints(Vector* point)
{
    std::vector<double> rbf;
    for (int i = 0; i < d_parameter_points.size(); i++)
    {
        rbf.push_back(obtainRBF(point, d_parameter_points[i]));
    }
    return rbf;
}

double Interpolator::obtainRBF(Vector* point1, Vector* point2)
{
    Vector diff;
    point1->minus(*point2, diff);
    double eps_norm_squared = d_epsilon * d_epsilon * diff.norm2();
    double res = 0.0;

    // Gaussian RBF
    if (d_rbf == "G")
    {
        res = std::exp(-eps_norm_squared);
    }
    // Multiquadric RBF
    else if (d_rbf == "MQ")
    {
        res = std::sqrt(1.0 + eps_norm_squared);
    }
    // Inverse quadratic RBF
    else if (d_rbf == "IQ")
    {
        res = 1.0 / (1.0 + eps_norm_squared);
    }
    // Inverse multiquadric RBF
    else if (d_rbf == "IMQ")
    {
        res = 1.0 / std::sqrt(1.0 + eps_norm_squared);
    }

    return res;
}

std::vector<Matrix*> obtainRotationMatrices(std::vector<Vector*> parameter_points,
        std::vector<Matrix*> bases,
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

    std::vector<Matrix*> rotation_matrices;

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
            Matrix* identity_matrix = new Matrix(bases[i]->numColumns(), bases[i]->numColumns(), false);
            for (int j = 0; j < identity_matrix->numColumns(); j++) {
                identity_matrix->item(j, j) = 1.0;
            }
            rotation_matrices.push_back(identity_matrix);
            continue;
        }

        Matrix* basis_mult_basis = bases[i]->transposeMult(bases[ref_point]);
        Matrix* basis = new Matrix(basis_mult_basis->numRows(), basis_mult_basis->numColumns(), false);
        Matrix* basis_right = new Matrix(basis_mult_basis->numColumns(), basis_mult_basis->numColumns(), false);

        // We need to compute the SVD of basis_mult_basis. Since it is
        // undistributed due to the transposeMult, let's use lapack's serial SVD
        // on rank 0 only.
        if (rank == 0)
        {
            Vector* sv = new Vector(basis_mult_basis->numColumns(), false);
            SerialSVD(basis_mult_basis, basis, sv, basis_right);
            delete sv;
        }

        // Broadcast U and V which are computed only on root.
        MPI_Bcast(basis->getData(), basis->numRows() * basis->numColumns(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(basis_right->getData(), basis_right->numRows() * basis_right->numColumns(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Obtain the rotation matrix.
        Matrix* rotation_matrix = basis->mult(basis_right);

        delete basis_mult_basis;
        delete basis;
        delete basis_right;
        rotation_matrices.push_back(rotation_matrix);
    }

    return rotation_matrices;
}

}
