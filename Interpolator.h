/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the Interpolator algorithm on the given snapshot matrix.

#ifndef included_Interpolator_h
#define included_Interpolator_h

#include <vector>

namespace CAROM {

class Matrix;
class Vector;

/**
 * Interpolator is a uninstantiable protected class that retains common
 * functionality between the MatrixInterpolator and VectorInterpolator classes.
 * The interpolation algorithm was adapted from "Gradient-based
 * Constrained Optimization Using a Database of Linear Reduced-Order Models"
 * by Y. Choi et al.
 */
class Interpolator
{

protected:

    /**
     * @brief Constructor.
     *
     * @param[in] parameter_points The parameter points.
     * @param[in] rotation_matrices The rotation matrices associated with
     *                              each parameter point.
     * @param[in] ref_point The index within the vector of parameter points
     *                      to the reference point
     */
    Interpolator(std::vector<Vector*> parameter_points,
                 std::vector<Matrix*> rotation_matrices,
                 int ref_point,
                 double epsilon = 1.0);

    /**
     * @brief The rank of the process this object belongs to.
     */
    int d_rank;

    /**
     * @brief The number of processors being run on.
     */
    int d_num_procs;

    /**
     * @brief The index within the vector of parameter points
     *        to the reference point
     */
    int d_ref_point;

    /**
     * @brief The RBF parameter that determines the width of influence.
     *        a small epsilon: larger influential width
     *        a large epsilon: smaller influential width
     */
    int d_epsilon;

    /**
     * @brief The sampled parameter points.
     */
    std::vector<Vector*> d_parameter_points;

    /**
     * @brief The reduced bases with compatible coordinates.
     */
    std::vector<Matrix*> d_rotation_matrices;

    /**
     * @brief The reduced elements in tangential space.
     */
    Matrix* d_lambda_T;

    /**
     * @brief Compute the RBF from the parameter points with the
     *        unsampled parameter point.
     *
     * @param[in] point The unsampled parameter point.
     */
    std::vector<double> obtainRBF(Vector* point);

private:

    /**
     * @brief Unimplemented default constructor.
     */
    Interpolator();

    /**
     * @brief Unimplemented copy constructor.
     */
    Interpolator(
        const Interpolator& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    Interpolator&
    operator = (
        const Interpolator& rhs);
};

/**
 * @brief Obtain the rotation matrices for all the parameter points using
 *        the basis of the reference point.
 *
 * @param[in] parameter_points The parameter points.
 * @param[in] bases The spatial basis associated with
 *                  each parameter point.
 * @param[in] ref_point The index within the vector of parameter points
 *                      to the reference point
 */
std::vector<Matrix*> obtainRotationMatrices(std::vector<Vector*> parameter_points,
        std::vector<Matrix*> bases,
        int ref_point);

}

#endif
