/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
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
#include <string>

namespace CAROM {

class Matrix;
class Vector;

/**
 * Interpolator is an uninstantiable protected class that retains common
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
     * @param[in] parameter_points  The parameter points.
     * @param[in] rotation_matrices The rotation matrices associated with
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
    Interpolator(std::vector<Vector*> parameter_points,
                 std::vector<Matrix*> rotation_matrices,
                 int ref_point,
                 std::string rbf,
                 std::string interp_method,
                 double closest_rbf_val = 0.9);

    ~Interpolator();

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
     * @brief The RBF type (gaussian, multiquadric, inverse quadratic, inverse
     *        multiquadric)
     */
    std::string d_rbf;

    /**
     * @brief The interpolation method (linear solve, inverse distance weighting,
     *        lagrangian polynomials)
     */
    std::string d_interp_method;

    /**
     * @brief The RBF parameter that determines the width of influence.
     *        a small epsilon: larger influential width
     *        a large epsilon: smaller influential width
     */
    double d_epsilon;

    /**
     * @brief The sampled parameter points.
     */
    std::vector<Vector*> d_parameter_points;

    /**
     * @brief The reduced bases with compatible coordinates.
     */
    std::vector<Matrix*> d_rotation_matrices;

    /**
     * @brief The RHS of the linear solve in tangential space.
     */
    Matrix* d_lambda_T;

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
 * @brief Compute the RBF from the parameter points with the
 *        unsampled parameter point.
 *
 * @param[in] parameter_points The parameter points.
 * @param[in] interp_method  The interpolation method type
 *                           ("LS" == linear solve,
 *                           "IDW" == inverse distance weighting,
 *                           "LP" == lagrangian polynomials)
 * @param[in] rbf Which RBF to compute.
 * @param[in] epsilon   The RBF parameter that determines the width of
                        influence.
 * @param[in] point The unsampled parameter point.
 */
std::vector<double> obtainRBFToTrainingPoints(std::vector<Vector*>
        parameter_points,
        std::string interp_method, std::string rbf, double epsilon, Vector* point);

/**
 * @brief Compute the sum of the RBF weights.
 *
 * @param[in] rbf The vector holding the rbfs of the training points.
 */
double rbfWeightedSum(std::vector<double>& rbf);

/**
 * @brief Compute the RBF between two points.
 *
 * @param[in] rbf Which RBF to compute.
 * @param[in] epsilon   The RBF parameter that determines the width of
                        influence.
 * @param[in] point1 The first point.
 * @param[in] point2 The second point.
 */
double obtainRBF(std::string rbf, double epsilon, Vector* point1,
                 Vector* point2);

/**
 * @brief Convert closest RBF value to an epsilon value.
 *
 * @param[in] parameter_points  The parameter points.
 * @param[in] rbf               Which RBF to compute.
 * @param[in] closest_rbf_val   The RBF parameter determines the width of influence.
 *                              Set the RBF value of the nearest two parameter points to a value between 0.0 to 1.0
 */
double convertClosestRBFToEpsilon(std::vector<Vector*> parameter_points,
                                  std::string rbf, double closest_rbf_val);

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
std::vector<Matrix*> obtainRotationMatrices(std::vector<Vector*>
        parameter_points,
        std::vector<Matrix*> bases,
        int ref_point);

}

#endif
