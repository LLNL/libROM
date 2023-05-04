/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the VectorInterpolator algorithm on the given snapshot matrix.

#ifndef included_VectorInterpolator_h
#define included_VectorInterpolator_h

#include "Interpolator.h"
#include <vector>
#include <string>

namespace CAROM {

class Matrix;
class Vector;

/**
 * VectorInterpolator interpolates reduced vectors of a set of parameter points
 * and returns an interpolated reduced vector for an unseen parameter point.
 * The interpolation algorithm was adapted from "Gradient-based
 * Constrained Optimization Using a Database of Linear Reduced-Order Models"
 * by Y. Choi et al.
 */
class VectorInterpolator : public Interpolator
{
public:

    /**
     * @brief Constructor.
     *
     * @param[in] parameter_points  The parameter points.
     * @param[in] rotation_matrices The rotation matrices associated with
     *                              each parameter point.
     * @param[in] reduced_vectors   The reduced vectors associated with
     *                              each parameter point.
     * @param[in] ref_point         The index within the vector of parameter points
     *                              to the reference point
     * @param[in] rbf               The RBF type ("G" == gaussian,
     *                              "IQ" == inverse quadratic,
     *                              "IMQ" == inverse multiquadric)
     * @param[in] interp_method     The interpolation method type
     *                              ("LS" == linear solve,
     *                              "IDW" == inverse distance weighting,
     *                              "LP" == lagrangian polynomials)
     * @param[in] closest_rbf_val   The RBF parameter determines the width of influence.
     *                              Set the RBF value of the nearest two parameter points to a value between 0.0 to 1.0
     */
    VectorInterpolator(std::vector<Vector*> parameter_points,
                       std::vector<Matrix*> rotation_matrices,
                       std::vector<Vector*> reduced_vectors,
                       int ref_point,
                       std::string rbf = "G",
                       std::string interp_method = "LS",
                       double closest_rbf_val = 0.9);

    ~VectorInterpolator();

    /**
     * @brief Obtain the interpolated reduced vector of the unsampled parameter point.
     *
     * @param[in] point The unsampled parameter point.
     */
    Vector* interpolate(Vector* point);

private:

    /**
     * @brief Unimplemented default constructor.
     */
    VectorInterpolator();

    /**
     * @brief Unimplemented copy constructor.
     */
    VectorInterpolator(
        const VectorInterpolator& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    VectorInterpolator&
    operator = (
        const VectorInterpolator& rhs);

    /**
     * @brief Solve the system of equations of the gammas to obtain the
     *        lambda for the P matrix.
     */
    void obtainLambda();

    /**
     * @brief Obtain the interpolated vector of the unsampled parameter point
     *           in log space.
     *
     * @param[in] rbf The RBF values between the parameter points and
     *                the unsampled parameter point.
     */
    Vector* obtainLogInterpolatedVector(std::vector<double>& rbf);

    /**
     * @brief The reduced vectors with compatible coordinates.
     */
    std::vector<Vector*> d_rotated_reduced_vectors;

    /**
     * @brief The reduced elements in tangential space.
     */
    std::vector<Vector*> d_gammas;
};

Vector* obtainInterpolatedVector(std::vector<Vector*> data, Matrix* f_T,
                                 std::string interp_method, std::vector<double>& rbf);

Matrix* solveLinearSystem(std::vector<Vector*> parameter_points,
                          std::vector<Vector*> data, std::string interp_method,
                          std::string rbf, double& epsilon);

}

#endif
