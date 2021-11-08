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
        SLPK_Matrix svd_input;

        int *row_offset = new int[num_procs + 1];
        row_offset[num_procs] = basis_mult_basis->numDistributedRows();
        row_offset[rank] = basis_mult_basis->numRows();

        CAROM_VERIFY(MPI_Allgather(MPI_IN_PLACE,
                                   1,
                                   MPI_INT,
                                   row_offset,
                                   1,
                                   MPI_INT,
                                   MPI_COMM_WORLD) == MPI_SUCCESS);
        for (int j = num_procs - 1; j >= 0; j--) {
            row_offset[j] = row_offset[j + 1] - row_offset[j];
        }

        CAROM_VERIFY(row_offset[0] == 0);

        int d_blocksize = row_offset[num_procs] / num_procs;
        if (row_offset[num_procs] % num_procs != 0) d_blocksize += 1;

        initialize_matrix(&svd_input, basis_mult_basis->numColumns(),
                          basis_mult_basis->numDistributedRows(),
                          1, num_procs, d_blocksize, d_blocksize);

        for (int rank = 0; rank < num_procs; ++rank)
        {
            scatter_block(&svd_input, 1, row_offset[rank] + 1,
                          basis_mult_basis->getData(),
                          basis_mult_basis->numColumns(),
                          row_offset[rank + 1] - row_offset[rank], rank);
        }

        std::unique_ptr<SVDManager> d_factorizer(new SVDManager);

        // This block does the actual ScaLAPACK call to do the factorization.
        svd_init(d_factorizer.get(), &svd_input);
        d_factorizer->dov = 1;
        factorize(d_factorizer.get());
        free_matrix_data(&svd_input);

        // Allocate the appropriate matrices and gather their elements.
        Matrix* basis = new Matrix(basis_mult_basis->numRows(), basis_mult_basis->numColumns(), basis_mult_basis->distributed());
        Matrix* basis_right = new Matrix(basis_mult_basis->numColumns(), basis_mult_basis->numColumns(), basis_mult_basis->distributed());
        for (int rank = 0; rank < num_procs; ++rank) {
            // V is computed in the transposed order so no reordering necessary.
            gather_block(&basis->item(0, 0), d_factorizer->V,
                         1, row_offset[static_cast<unsigned>(rank)]+1,
                         basis_mult_basis->numColumns(), row_offset[static_cast<unsigned>(rank) + 1] -
                         row_offset[static_cast<unsigned>(rank)],
                         rank);

            // U is computed in the transposed order so no reordering necessary.
            gather_block(&basis_right->item(0, 0), d_factorizer->U, 1, 1,
                         basis_mult_basis->numColumns(), basis_mult_basis->numColumns(), rank);
        }

        Matrix* rotation_matrix = basis->mult(basis_right);

        delete [] row_offset;
        delete basis_mult_basis;
        delete basis;
        delete basis_right;
        rotation_matrices.push_back(rotation_matrix);
    }

    return rotation_matrices;
}

}
