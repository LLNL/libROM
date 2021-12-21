/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the AdaptiveDMD algorithm on the given snapshot matrix. The
//              implemented dynamic mode decomposition algorithm is derived from
//              Tu et. al's paper "On Dynamic Mode Decomposition: Theory and
//              Applications": https://arxiv.org/abs/1312.0041

#ifndef included_AdaptiveDMD_h
#define included_AdaptiveDMD_h

#include "DMD.h"
#include <vector>
#include <complex>

namespace CAROM {

class Matrix;
class Vector;

/**
 * Class AdaptiveDMD implements the AdaptiveDMD algorithm on a given snapshot matrix.
 */
class AdaptiveDMD : public DMD
{
public:

    /**
     * @brief Constructor.
     *
     * @param[in] dim The full-order state dimension.
     * @param[in] adaptive_dt Whether the dt is adaptive.
     */
    AdaptiveDMD(int dim);

    /**
     * @brief Sample the new state, u_in.
     *
     * @pre u_in != 0
     * @pre t >= 0.0
     *
     * @param[in] u_in The new state.
     * @param[in] t    The time of the newly sampled state.
     */
    void takeSample(double* u_in, double t);

    /**
     * @brief Predict state given a new initial condition and time.
     *        The initial condition must be projected using projectInitialCondition
     *        for correct results.
     *
     * @param[in] init The initial condition.
     * @param[in] n The time of the outputted state (t/dt)
     */
    Vector* predict(const std::pair<Vector*, Vector*> init, double n);

private:

    /**
     * @brief Unimplemented default constructor.
     */
    AdaptiveDMD();

    /**
     * @brief Unimplemented copy constructor.
     */
    AdaptiveDMD(
        const AdaptiveDMD& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    AdaptiveDMD&
    operator = (
        const AdaptiveDMD& rhs);

    /**
     * @brief The stored times of each sample.
     */
    std::vector<double> d_sampled_times;

    /**
     * @brief Interpolate n to obtain the sampled time to input into the DMD
     *        prediction algorithm.
     *
     * @param[in] n The time of the outputted state (t/dt)
     */
    double interpolateSampledTime(double n);
};

}

#endif
