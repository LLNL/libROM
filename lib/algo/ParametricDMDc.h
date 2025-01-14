/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the ParametricDMDc algorithm to obtain DMDc model interpolant
//              at desired parameter point by interpolation of DMDc models at training parameter points.
//              The implemented dynamic mode decomposition algorithm is derived from
//              Tu et. al's paper "On Dynamic Mode Decomposition: Theory and
//              Applications": https://arxiv.org/abs/1312.0041
//              The interpolation algorithm was adapted from "Gradient-based
//              Constrained Optimization Using a Database of Linear Reduced-Order Models"
//              by Y. Choi et al. We extend this algorithm to DMDc.

#ifndef included_ParametricDMDc_h
#define included_ParametricDMDc_h

#include "manifold_interp/MatrixInterpolator.h"
#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include "mpi.h"

#include <complex>

namespace CAROM {

/**
 * @brief Constructor.
 *
 * @param[in] parametric_dmdc       The interpolant DMDc model at the desired point.
 * @param[in] parameter_points      The training parameter points.
 * @param[in] dmdcs                 The DMDc objects associated with
 *                                  each training parameter point.
 * @param[in] controls              The matrices of controls from previous runs
 *                                  which we use to interpolate.
 * @param[in] controls_interpolated The interpolated controls.
 * @param[in] desired_point         The desired point at which to create a parametric DMDc.
 * @param[in] rbf                   The RBF type ("G" == gaussian,
 *                                  "IQ" == inverse quadratic,
 *                                  "IMQ" == inverse multiquadric)
 * @param[in] interp_method         The interpolation method type
 *                                  ("LS" == linear solve,
 *                                  "IDW" == inverse distance weighting,
 *                                  "LP" == lagrangian polynomials)
 * @param[in] closest_rbf_val       The RBF parameter determines the width of influence.
 *                                  Set the RBF value of the nearest two parameter points
 *                                  to a value between 0.0 to 1.0
 * @param[in] reorthogonalize_W     Whether to reorthogonalize the interpolated W (basis) matrix.
 */
template <class T>
void getParametricDMDc(std::unique_ptr<T>& parametric_dmdc,
                       const std::vector<Vector>& parameter_points,
                       std::vector<T*>& dmdcs,
                       std::vector<std::shared_ptr<Matrix>> & controls,
                       std::shared_ptr<Matrix> & controls_interpolated,
                       const Vector & desired_point,
                       std::string rbf = "G",
                       std::string interp_method = "LS",
                       double closest_rbf_val = 0.9,
                       bool reorthogonalize_W = false)
{
    CAROM_VERIFY(parameter_points.size() == dmdcs.size());
    CAROM_VERIFY(dmdcs.size() > 1);
    for (int i = 0; i < dmdcs.size() - 1; i++)
    {
        CAROM_VERIFY(dmdcs[i]->d_dt == dmdcs[i + 1]->d_dt);
        CAROM_VERIFY(dmdcs[i]->d_k == dmdcs[i + 1]->d_k);
    }
    CAROM_VERIFY(closest_rbf_val >= 0.0 && closest_rbf_val <= 1.0);

    int mpi_init, rank;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<std::shared_ptr<CAROM::Matrix>> bases;
    std::vector<std::shared_ptr<CAROM::Matrix>> A_tildes;
    std::vector<std::shared_ptr<CAROM::Matrix>> B_tildes;

    for (int i = 0; i < dmdcs.size(); i++)
    {
        bases.push_back(dmdcs[i]->d_basis);
        A_tildes.push_back(dmdcs[i]->d_A_tilde);
        B_tildes.push_back(dmdcs[i]->d_B_tilde);
    }

    int ref_point = getClosestPoint(parameter_points, desired_point);
    std::vector<std::shared_ptr<CAROM::Matrix>> rotation_matrices =
                obtainRotationMatrices(parameter_points, bases, ref_point);

    CAROM::MatrixInterpolator basis_interpolator(parameter_points,
            rotation_matrices, bases, ref_point, "B", rbf, interp_method, closest_rbf_val);
    std::shared_ptr<CAROM::Matrix> W = basis_interpolator.interpolate(desired_point,
                                       reorthogonalize_W);

    CAROM::MatrixInterpolator A_tilde_interpolator(parameter_points,
            rotation_matrices, A_tildes, ref_point, "R", rbf, interp_method,
            closest_rbf_val);
    std::shared_ptr<CAROM::Matrix> A_tilde = A_tilde_interpolator.interpolate(
                desired_point);

    CAROM::MatrixInterpolator B_tilde_interpolator(parameter_points,
            rotation_matrices, B_tildes, ref_point, "R", rbf, interp_method,
            closest_rbf_val);
    std::shared_ptr<CAROM::Matrix> B_tilde = B_tilde_interpolator.interpolate(
                desired_point);

    CAROM::MatrixInterpolator control_interpolator(parameter_points,
            rotation_matrices, controls, ref_point, "R", rbf, interp_method,
            closest_rbf_val);

    controls_interpolated = control_interpolator.interpolate(desired_point);

    // Calculate the right eigenvalues/eigenvectors of A_tilde
    ComplexEigenPair eigenpair = NonSymmetricRightEigenSolve(*A_tilde);
    std::vector<std::complex<double>> eigs = eigenpair.eigs;

    // Calculate phi (phi = W * eigenvectors)
    std::shared_ptr<CAROM::Matrix> phi_real = W->mult(*eigenpair.ev_real);
    std::shared_ptr<CAROM::Matrix> phi_imaginary = W->mult(*eigenpair.ev_imaginary);

    parametric_dmdc.reset(new T(eigs, phi_real, phi_imaginary, B_tilde,
                                dmdcs[0]->d_k,dmdcs[0]->d_dt,
                                dmdcs[0]->d_t_offset, dmdcs[0]->d_state_offset, W));

    delete eigenpair.ev_real;
    delete eigenpair.ev_imaginary;
}

/**
 * @brief Constructor.
 *
 * @param[in] parameter_points      The parameter points.
 * @param[in] dmdc_paths            The paths to the saved DMD objects associated with
 *                                  each parameter point.
 * @param[in] desired_point         The desired point at which to create a parametric DMD.
 * @param[in] controls              The matrices of controls from previous runs which we
 *                                  use to interpolate.
 * @param[in] controls_interpolated The interpolated controls.
 * @param[in] rbf                   The RBF type ("G" == gaussian,
 *                                  "IQ" == inverse quadratic,
 *                                  "IMQ" == inverse multiquadric)
 * @param[in] interp_method         The interpolation method type
 *                                  ("LS" == linear solve,
 *                                  "IDW" == inverse distance weighting,
 *                                  "LP" == lagrangian polynomials)
 * @param[in] closest_rbf_val       The RBF parameter determines the width of influence.
 *                                  Set the RBF value of the nearest two parameter points
 *                                  to a value between 0.0 to 1.0
 * @param[in] reorthogonalize_W     Whether to reorthogonalize the interpolated W (basis) matrix.
 */
template <class T>
void getParametricDMDc(std::unique_ptr<T>& parametric_dmdc,
                       const std::vector<Vector>& parameter_points,
                       std::vector<std::string>& dmdc_paths,
                       std::vector<std::shared_ptr<Matrix>> & controls,
                       std::shared_ptr<Matrix> & controls_interpolated,
                       const Vector & desired_point,
                       std::string rbf = "G",
                       std::string interp_method = "LS",
                       double closest_rbf_val = 0.9,
                       bool reorthogonalize_W = false)
{
    std::vector<T*> dmdcs;
    for (int i = 0; i < dmdc_paths.size(); i++)
    {
        T* dmdc = new T(dmdc_paths[i]);
        dmdcs.push_back(dmdc);
    }

    getParametricDMDc(parametric_dmdc, parameter_points, dmdcs, controls,
                      controls_interpolated,
                      desired_point,
                      rbf, interp_method, closest_rbf_val,
                      reorthogonalize_W);

    for (int i = 0; i < dmdcs.size(); i++)
    {
        delete dmdcs[i];
    }
}

}

#endif
