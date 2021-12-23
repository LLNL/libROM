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
//              between the samples. This algorithm interpolates the ununiformly
//              sampled samples such that they are uniformly distanced and the
//              dt is constant (a prequisite of the DMD algorithm). The smaller
//              dt is, the finer the fidelity of the interpolation.

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
     * @param[in] desired_dt  The dt to uniformly interpolate the samples between.
     * @param[in] interp_method  The interpolation method type ("LS" == linear solve,
     *                           "IDW" == inverse distance weighting, "LP" == lagrangian polynomials)
     * @param[in] rbf         Which RBF to compute.
     */
    AdaptiveDMD(int dim, double desired_dt, std::string interp_method, std::string rbf);

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
     * @brief Internal function to obtain the interpolated snapshots.
     */
    const Matrix* interpolateSnapshots();

    /**
     * @brief The stored times of each sample.
     */
    std::vector<Vector*> d_sampled_times;

    /**
     * @brief The dt to uniformly interpolate the samples between.
     */
    double d_dt;

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
};

}

#endif
