/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the MatrixInterpolater algorithm on the given snapshot matrix.

#ifndef included_MatrixInterpolater_h
#define included_MatrixInterpolater_h

#include "Interpolater.h"
#include <vector>
#include <string>

namespace CAROM {

class Matrix;
class Vector;

class MatrixInterpolater : public Interpolater
{
public:

    /**
     * @brief Constructor.
     *
     * @param[in] parameter_points The parameter points.
     * @param[in] rotation_matrices The rotation matrices associated with
     *                              each parameter point.
     * @param[in] reduced_matrices The reduced matrices associated with
     *                             each parameter point.
     * @param[in] ref_point The index within the vector of parameter points
     *                      to the reference point
     * @param[in] matrix_type The type of matrix (R = real, B = basis,
     *                        NS = nonsingular, SPD = symmetric positive-definite)
     */
    MatrixInterpolater(std::vector<Vector*> parameter_points,
                       std::vector<Matrix*> rotation_matrices,
                       std::vector<Matrix*> reduced_matrices,
                       int ref_point,
                       std::string matrix_type);

    /**
     * @brief Obtain the interpolated reduced matrix of the unsampled parameter point.
     *
     * @param[in] point The unsampled parameter point.
     */
    Matrix* interpolate(Vector* point);

private:

    /**
     * @brief Unimplemented default constructor.
     */
    MatrixInterpolater();

    /**
     * @brief Unimplemented copy constructor.
     */
    MatrixInterpolater(
        const MatrixInterpolater& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    MatrixInterpolater&
    operator = (
        const MatrixInterpolater& rhs);

    /**
     * @brief Solve the system of equations of the gammas to obtain the
     *        lambda for the P matrix.
     *
     * @param[in] gammas The vector of gamma matrices.
     */
    void obtainLambda(std::vector<Matrix*> gammas);

    /**
     * @brief Obtain the interpolated matrix of the unsampled parameter point
     *        in log space.
     * @param[in] inv_q The inverse quadratic RBF values between the
     *                  parameter points and the unsampled parameter point.
     */
    Matrix* obtainLogInterpolatedMatrix(std::vector<double> inv_q);

    /**
     * @brief Obtain the interpolated SPD reduced matrix for the unsampled
     *        parameter point.
     *
     * @param[in] point The unsampled parameter point.
     */
    Matrix* interpolateSPDMatrix(Vector* point);

    /**
     * @brief Obtain the interpolated non-singular reduced matrix for the
     *        unsampled parameter point.
     *
     * @param[in] point The unsampled parameter point.
     */
    Matrix* interpolateNonSingularMatrix(Vector* point);

    /**
     * @brief Obtain the interpolated reduced matrix for the
     *        unsampled parameter point.
     *
     * @param[in] point The unsampled parameter point.
     */
    Matrix* interpolateMatrix(Vector* point);

    /**
     * @brief The matrix type/manifold
     */
    std::string d_matrix_type;

    /**
     * @brief The reduced matrices with compatible coordinates.
     */
    std::vector<Matrix*> d_rotated_reduced_matrices;

    /**
     * @brief The reduced matrix of the reference point to the half power.
     */
    Matrix* d_x_half_power;
};

}

#endif
