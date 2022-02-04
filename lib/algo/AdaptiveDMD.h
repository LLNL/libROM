/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the AdaptiveDMD algorithm on the given snapshot matrix.
//              The AdaptiveDMD algorithm should be used if the dt changes
//              between the samples. This algorithm uniformly interpolates the samples
//              that may have been taken with variable steps, using the constant step dt
//              (a prequisite of the DMD algorithm). The smaller dt is, the finer
//              the fidelity of the interpolation. This algorithm also works in
//              the case that the first sample does not start from t = 0.0 by
//              incorporating a time offset.

#ifndef included_AdaptiveDMD_h
#define included_AdaptiveDMD_h

#include "DMD.h"
#include <vector>

namespace CAROM {

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
     * @param[in] dim         The full-order state dimension.
     * @param[in] desired_dt  The constant step size for uniform interpolation of samples.
     * @param[in] rbf         The RBF type ("G" == gaussian, "MQ" == multiquadric,
     *                        "IQ" == inverse quadratic, "IMQ" == inverse
     *                        multiquadric)
     * @param[in] interp_method  The interpolation method type ("LS" == linear solve,
     *                           "IDW" == inverse distance weighting, "LP" == lagrangian polynomials)
     * @param[in] epsilon   The RBF parameter that determines the width of
                            influence.
     */
    AdaptiveDMD(int dim, double desired_dt, std::string rbf = "G", std::string interp_method = "LS", double epsilon = -1.0);

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
     * @param[in] energy_fraction The energy fraction to keep after doing SVD.
     */
    void train(double energy_fraction);

    /**
     * @param[in] k The number of modes (eigenvalues) to keep after doing SVD.
     */
    void train(int k);

    /**
     * @brief Get the true dt between interpolated snapshots. 
     */
    double getTrueDt();

    /**
     * @brief Get the interpolated snapshot matrix contained within d_interp_snapshots.
     */
    const Matrix* getInterpolatedSnapshots(); 

    /**
     * @brief Predict state given a time. Uses the projected initial condition of the
     *        training dataset (the first column).
     *
     * @param[in] t The time of the outputted state
     */
    Vector* predict(double t);

    /**
     * @brief Predict state given a new initial condition and time.
     *        The initial condition must be projected using projectInitialCondition
     *        for correct results.
     *
     * @param[in] init The initial condition.
     * @param[in] t The time of the outputted state
     */
    Vector* predict(const std::pair<Vector*, Vector*> init, double t);

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
    std::vector<Vector*> d_sampled_times;

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
     * @brief std::vector holding the interpolated snapshots.
     */
    std::vector<Vector*> d_interp_snapshots;

    /**
     * @brief Internal function to obtain the interpolated snapshots.
     */
    void interpolateSnapshots();

    /**
     * @brief The RBF parameter that determines the width of influence.
     *        a small epsilon: larger influential width
     *        a large epsilon: smaller influential width
     */
    double d_epsilon;

    /**
     * @brief The time offset of the first sample.
     */
    double d_t_offset;
};

}

#endif
