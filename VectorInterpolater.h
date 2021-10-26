/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the VectorInterpolater algorithm on the given snapshot matrix.

#ifndef included_VectorInterpolater_h
#define included_VectorInterpolater_h

#include "Interpolater.h"
#include <vector>
#include <string>

namespace CAROM {

class Matrix;
class Vector;

/**
 * VectorInterpolater interpolates reduced vectors of a set of parameter points
 * and returns an interpolated reduced vector for an unseen parameter point.
 */
class VectorInterpolater : public Interpolater
{
public:

    /**
     * @brief Constructor.
     *
     * @param[in] parameter_points The parameter points.
     * @param[in] rotation_matrices The rotation matrices associated with
     *                              each parameter point.
     * @param[in] reduced_vectors The reduced vectors associated with
     *                            each parameter point.
     * @param[in] ref_point The index within the vector of parameter points
     *                      to the reference point
     */
    VectorInterpolater(std::vector<Vector*> parameter_points,
                       std::vector<Matrix*> rotation_matrices,
                       std::vector<Vector*> reduced_vectors,
                       int ref_point);

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
    VectorInterpolater();

    /**
     * @brief Unimplemented copy constructor.
     */
    VectorInterpolater(
        const VectorInterpolater& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    VectorInterpolater&
    operator = (
        const VectorInterpolater& rhs);

    /**
     * @brief Solve the system of equations of the gammas to obtain the
     *        lambda for the P matrix.
     *
     * @param[in] gammas The vector of gamma matrices.
     */
    void obtainLambda(std::vector<Vector*> gammas);

    /**
     * @brief Obtain the interpolated vector of the unsampled parameter point
     *           in log space.
     *
     * @param[in] inv_q The inverse quadratic RBF values between the
     *                  parameter points and the unsampled parameter point.
     */
    Vector* obtainLogInterpolatedVector(std::vector<double> inv_q);

    /**
     * @brief The reduced vectors with compatible coordinates.
     */
    std::vector<Vector*> d_rotated_reduced_vectors;
};

}

#endif
