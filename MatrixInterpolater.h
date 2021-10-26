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
     * @param[in] dim The full-order state dimension.
     */
    MatrixInterpolater(std::vector<Vector*> parameter_points,
                       std::vector<Matrix*> rotation_matrices,
                       std::vector<Matrix*> reduced_matrices,
                       int ref_point,
                       std::string matrix_type);

   /**
    * @brief Interpolate to unsampled parameter point.
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
         * @brief Interpolate to unsampled parameter point.
         *
         * @param[in] point The unsampled parameter point.
         */
        void obtainLambda(std::vector<Matrix*> gammas);

        /**
         * @brief Interpolate to unsampled parameter point.
         *
         * @param[in] point The unsampled parameter point.
         */
        Matrix* obtainLogInterpolatedMatrix(std::vector<double> inv_q);

        /**
         * @brief Interpolate to unsampled parameter point.
         *
         * @param[in] point The unsampled parameter point.
         */
         Matrix* interpolateSPDMatrix(Vector* point);

         /**
          * @brief Interpolate to unsampled parameter point.
          *
          * @param[in] point The unsampled parameter point.
          */
          Matrix* interpolateNonSingularMatrix(Vector* point);

          /**
           * @brief Interpolate to unsampled parameter point.
           *
           * @param[in] point The unsampled parameter point.
           */
           Matrix* interpolateMatrix(Vector* point);

    /**
     * @brief The matrix manifold
     */
    std::string d_matrix_type;

    /**
     * @brief The reduced matrices with compatible coordinates.
     */
    std::vector<Matrix*> d_rotated_reduced_matrices;

    Matrix* d_x_half_power;
};

}

#endif
