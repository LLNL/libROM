/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the Interpolater algorithm on the given snapshot matrix.

#ifndef included_Interpolater_h
#define included_Interpolater_h

#include <vector>

namespace CAROM {

class Matrix;
class Vector;

class Interpolater
{
public:

    /**
     * @brief Constructor.
     *
     * @param[in] dim The full-order state dimension.
     */
    Interpolater(std::vector<Vector*> parameter_points,
                 std::vector<Matrix*> rotation_matrices,
                 int ref_point);

protected:

    /**
     * @brief The rank of the process this object belongs to.
     */
    int d_rank;

    /**
     * @brief The number of processors being run on.
     */
    int d_num_procs;

    /**
     * @brief The reference point;
     */
    int d_ref_point;

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
     * @brief Interpolate to unsampled parameter point.
     *
     * @param[in] point The unsampled parameter point.
     */
    std::vector<double> obtainRBF(Vector* point);

private:

    /**
     * @brief Unimplemented default constructor.
     */
    Interpolater();

    /**
     * @brief Unimplemented copy constructor.
     */
    Interpolater(
        const Interpolater& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    Interpolater&
    operator = (
        const Interpolater& rhs);
};

std::vector<Matrix*> obtainRotationMatrices(std::vector<Vector*> parameter_points,
                                            std::vector<Matrix*> bases,
                                            int ref_point);

}

#endif
