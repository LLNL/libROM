/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Implementation of the ParametricDMD algorithm.

#include "ParametricDMD.h"
#include "manifold_interp/MatrixInterpolator.h"
#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include "mpi.h"

#include <complex>

namespace CAROM {

DMD* getParametricDMD(std::vector<Vector*>& parameter_points,
                      std::vector<DMD*>& dmds,
                      Vector* desired_point,
                      std::string rbf,
                      std::string interp_method,
                      double closest_rbf_val,
                      bool reorthogonalize_W)
{
    CAROM_VERIFY(parameter_points.size() == dmds.size());
    CAROM_VERIFY(dmds.size() > 1);
    for (int i = 0; i < dmds.size() - 1; i++)
    {
        CAROM_VERIFY(dmds[i]->d_dt == dmds[i + 1]->d_dt);
        CAROM_VERIFY(dmds[i]->d_t_offset == dmds[i + 1]->d_t_offset);
        CAROM_VERIFY(dmds[i]->d_k == dmds[i + 1]->d_k);
    }
    CAROM_VERIFY(closest_rbf_val >= 0.0 && closest_rbf_val <= 1.0);

    int mpi_init, rank;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<Matrix*> bases;
    std::vector<Matrix*> A_tildes;
    for (int i = 0; i < dmds.size(); i++)
    {
        bases.push_back(dmds[i]->d_basis);
        A_tildes.push_back(dmds[i]->d_A_tilde);
    }

    int ref_point = getClosestPoint(parameter_points, desired_point);
    std::vector<CAROM::Matrix*> rotation_matrices = obtainRotationMatrices(parameter_points,
            bases, ref_point);

    CAROM::MatrixInterpolator basis_interpolator(parameter_points,
            rotation_matrices, bases, ref_point, "B", rbf, interp_method, closest_rbf_val);
    Matrix* W = basis_interpolator.interpolate(desired_point, reorthogonalize_W);

    CAROM::MatrixInterpolator A_tilde_interpolator(parameter_points,
            rotation_matrices, A_tildes, ref_point, "R", rbf, interp_method, closest_rbf_val);
    Matrix* A_tilde = A_tilde_interpolator.interpolate(desired_point);

    // Calculate the right eigenvalues/eigenvectors of A_tilde
    ComplexEigenPair eigenpair = NonSymmetricRightEigenSolve(A_tilde);
    std::vector<std::complex<double>> eigs = eigenpair.eigs;

    // Calculate phi (phi = W * eigenvectors)
    Matrix* phi_real = W->mult(eigenpair.ev_real);
    Matrix* phi_imaginary = W->mult(eigenpair.ev_imaginary);

    DMD* desired_dmd = new DMD(eigs, phi_real, phi_imaginary, dmds[0]->d_k,
                               dmds[0]->d_dt, dmds[0]->d_t_offset);

    delete W;
    delete A_tilde;
    delete eigenpair.ev_real;
    delete eigenpair.ev_imaginary;

    return desired_dmd;
}

DMD* getParametricDMD(std::vector<Vector*>& parameter_points,
                      std::vector<std::string>& dmd_paths,
                      Vector* desired_point,
                      std::string rbf,
                      std::string interp_method,
                      double closest_rbf_val,
                      bool reorthogonalize_W)
{
    std::vector<DMD*> dmds;
    for (int i = 0; i < dmd_paths.size(); i++)
    {
        DMD* dmd = new DMD(dmd_paths[i]);
        dmds.push_back(dmd);
    }

    DMD* desired_dmd = getParametricDMD(parameter_points, dmds, desired_point,
                                        rbf, interp_method, closest_rbf_val,
                                        reorthogonalize_W);
    for (int i = 0; i < dmds.size(); i++)
    {
        delete dmds[i];
    }

    return desired_dmd;
}

}
