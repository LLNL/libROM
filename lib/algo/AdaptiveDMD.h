/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
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
//              the fidelity of the interpolation.

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
     * @param[in] dim              The full-order state dimension.
     * @param[in] desired_dt       The constant step size for uniform interpolation of samples.
     *                             If set equal to or below 0.0, desired_dt will be set to the median of
     *                             the different dt's between the samples.
     * @param[in] rbf              The RBF type ("G" == gaussian,
     *                             "IQ" == inverse quadratic,
     *                             "IMQ" == inverse multiquadric)
     * @param[in] interp_method    The interpolation method type
     *                             ("LS" == linear solve,
     *                             "IDW" == inverse distance weighting,
     *                             "LP" == lagrangian polynomials)
     * @param[in] closest_rbf_val  The RBF parameter determines the width of influence.
     *                             Set the RBF value of the nearest two parameter points
     *                             to a value between 0.0 to 1.0.
     * @param[in] alt_output_basis Whether to use the alternative basis for
     *                             output, i.e. phi = U^(+)*V*Omega^(-1)*X.
     * @param[in] state_offset     The state offset.
     */
    AdaptiveDMD(int dim, double desired_dt = -1.0, std::string rbf = "G",
                std::string interp_method = "LS",
                double closest_rbf_val = 0.9,
                bool alt_output_basis = false,
                Vector* state_offset = NULL);

    /**
     * @brief Destroy the AdaptiveDMD object
     */
    ~AdaptiveDMD();

    /**
     * @param[in] energy_fraction The energy fraction to keep after doing SVD.
     * @param[in] W0              The initial basis to prepend to W.
     * @param[in] linearity_tol   The tolerance for determining whether a column
                                  of W is linearly independent with W0.
     */
    void train(double energy_fraction, const Matrix* W0 = NULL,
               double linearity_tol = 0.0);

    /**
     * @param[in] k The number of modes (eigenvalues) to keep after doing SVD.
     * @param[in] W0              The initial basis to prepend to W.
     * @param[in] linearity_tol   The tolerance for determining whether a column
                                  of W is linearly independent with W0.
     */
    void train(int k, const Matrix* W0 = NULL, double linearity_tol = 0.0);

    /**
     * @brief Get the true dt between interpolated snapshots.
     */
    double getTrueDt() const;

    /**
     * @brief Get the interpolated snapshot matrix contained within d_interp_snapshots.
     */
    const Matrix* getInterpolatedSnapshots();

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
     * @brief The RBF value of the nearest two parameter points
     */
    double d_closest_rbf_val;
};

}

#endif
