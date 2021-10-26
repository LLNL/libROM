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

class VectorInterpolater : public Interpolater
{
public:

    /**
     * @brief Constructor.
     *
     * @param[in] dim The full-order state dimension.
     */
    VectorInterpolater(std::vector<Vector*> parameter_points,
                       std::vector<Matrix*> rotation_matrices,
                       std::vector<Vector*> reduced_vectors,
                       int ref_point);

   /**
    * @brief Interpolate to unsampled parameter point.
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
         * @brief Interpolate to unsampled parameter point.
         *
         * @param[in] point The unsampled parameter point.
         */
        void obtainLambda(std::vector<Vector*> gammas);

        /**
         * @brief Interpolate to unsampled parameter point.
         *
         * @param[in] point The unsampled parameter point.
         */
        Vector* obtainLogInterpolatedVector(std::vector<double> inv_q);

    /**
     * @brief The reduced vectors with compatible coordinates.
     */
    std::vector<Vector*> d_rotated_reduced_vectors;
};

}

#endif
