/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the ParametricDMD algorithm to obtain DMD model interpolant
//              at desired parameter point by interpolation of DMD models at training parameter points.
//              The implemented dynamic mode decomposition algorithm is derived from
//              Tu et. al's paper "On Dynamic Mode Decomposition: Theory and
//              Applications": https://arxiv.org/abs/1312.0041
//              The interpolation algorithm was adapted from "Gradient-based
//              Constrained Optimization Using a Database of Linear Reduced-Order Models"
//              by Y. Choi et al.

#ifndef included_ParametricDMD_h
#define included_ParametricDMD_h

#include "manifold_interp/MatrixInterpolator.h"
#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include "mpi.h"

#include <complex>

namespace CAROM {

/**
 * @brief Constructor.
 *
 * @param[in] parametric_dmd    The interpolant DMD model at the desired point.
 * @param[in] parameter_points  The training parameter points.
 * @param[in] dmds              The DMD objects associated with
 *                              each training parameter point.
 * @param[in] desired_point     The desired point at which to create a parametric DMD.
 * @param[in] rbf               The RBF type ("G" == gaussian,
 *                              "IQ" == inverse quadratic,
 *                              "IMQ" == inverse multiquadric)
 * @param[in] interp_method     The interpolation method type
 *                              ("LS" == linear solve,
 *                              "IDW" == inverse distance weighting,
 *                              "LP" == lagrangian polynomials)
 * @param[in] closest_rbf_val   The RBF parameter determines the width of influence.
 *                              Set the RBF value of the nearest two parameter points to a value between 0.0 to 1.0
 * @param[in] reorthogonalize_W Whether to reorthogonalize the interpolated W (basis) matrix.
 */
template <class T>
void getParametricDMD(T*& parametric_dmd,
                      std::vector<Vector*>& parameter_points,
                      std::vector<T*>& dmds,
                      Vector* desired_point,
                      std::string rbf = "G",
                      std::string interp_method = "LS",
                      double closest_rbf_val = 0.9,
                      bool reorthogonalize_W = false)
{
    CAROM_VERIFY(parameter_points.size() == dmds.size());
    CAROM_VERIFY(dmds.size() > 1);
    for (int i = 0; i < dmds.size() - 1; i++)
    {
        CAROM_VERIFY(dmds[i]->d_dt == dmds[i + 1]->d_dt);
        CAROM_VERIFY(dmds[i]->d_k == dmds[i + 1]->d_k);
    }
    CAROM_VERIFY(closest_rbf_val >= 0.0 && closest_rbf_val <= 1.0);

    int mpi_init, rank;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<CAROM::Matrix*> bases;
    std::vector<CAROM::Matrix*> A_tildes;
    for (int i = 0; i < dmds.size(); i++)
    {
        bases.push_back(dmds[i]->d_basis);
        A_tildes.push_back(dmds[i]->d_A_tilde);
    }

    int ref_point = getClosestPoint(parameter_points, desired_point);
    std::vector<CAROM::Matrix*> rotation_matrices = obtainRotationMatrices(
                parameter_points,
                bases, ref_point);

    CAROM::MatrixInterpolator basis_interpolator(parameter_points,
            rotation_matrices, bases, ref_point, "B", rbf, interp_method, closest_rbf_val);
    CAROM::Matrix* W = basis_interpolator.interpolate(desired_point,
                       reorthogonalize_W);

    CAROM::MatrixInterpolator A_tilde_interpolator(parameter_points,
            rotation_matrices, A_tildes, ref_point, "R", rbf, interp_method,
            closest_rbf_val);
    CAROM::Matrix* A_tilde = A_tilde_interpolator.interpolate(desired_point);

    // Calculate the right eigenvalues/eigenvectors of A_tilde
    ComplexEigenPair eigenpair = NonSymmetricRightEigenSolve(A_tilde);
    std::vector<std::complex<double>> eigs = eigenpair.eigs;

    // Calculate phi (phi = W * eigenvectors)
    CAROM::Matrix* phi_real = W->mult(eigenpair.ev_real);
    CAROM::Matrix* phi_imaginary = W->mult(eigenpair.ev_imaginary);

    parametric_dmd = new T(eigs, phi_real, phi_imaginary, dmds[0]->d_k,
                           dmds[0]->d_dt, dmds[0]->d_t_offset,
                           dmds[0]->d_state_offset);

    delete W;
    delete A_tilde;
    delete eigenpair.ev_real;
    delete eigenpair.ev_imaginary;

    for (auto m : rotation_matrices)
        delete m;
}

/**
 * @brief Constructor.
 *
 * @param[in] parameter_points  The parameter points.
 * @param[in] dmd_paths         The paths to the saved DMD objects associated with
 *                              each parameter point.
 * @param[in] desired_point     The desired point at which to create a parametric DMD.
 * @param[in] rbf               The RBF type ("G" == gaussian,
 *                              "IQ" == inverse quadratic,
 *                              "IMQ" == inverse multiquadric)
 * @param[in] interp_method     The interpolation method type
 *                              ("LS" == linear solve,
 *                              "IDW" == inverse distance weighting,
 *                              "LP" == lagrangian polynomials)
 * @param[in] closest_rbf_val   The RBF parameter determines the width of influence.
 *                              Set the RBF value of the nearest two parameter points to a value between 0.0 to 1.0
 * @param[in] reorthogonalize_W Whether to reorthogonalize the interpolated W (basis) matrix.
 */
template <class T>
void getParametricDMD(T*& parametric_dmd,
                      std::vector<Vector*>& parameter_points,
                      std::vector<std::string>& dmd_paths,
                      Vector* desired_point,
                      std::string rbf = "G",
                      std::string interp_method = "LS",
                      double closest_rbf_val = 0.9,
                      bool reorthogonalize_W = false)
{
    std::vector<T*> dmds;
    for (int i = 0; i < dmd_paths.size(); i++)
    {
        T* dmd = new T(dmd_paths[i]);
        dmds.push_back(dmd);
    }

    getParametricDMD(parametric_dmd, parameter_points, dmds, desired_point,
                     rbf, interp_method, closest_rbf_val,
                     reorthogonalize_W);
    for (int i = 0; i < dmds.size(); i++)
    {
        delete dmds[i];
    }
}

}

#endif
