/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Computes the DMD algorithm on the given snapshot matrix 
//              with non-uniform sampling time steps. 
//              The implemented dynamic mode decomposition algorithm is modified from
//              Tu et. al's paper "On Dynamic Mode Decomposition: Theory and
//              Applications": https://arxiv.org/abs/1312.0041
//              Instead of approximating the discrete dynamics, this algorithm 
//              approximates the continuous dynamics linearly. 
//              This algorithm also works in the case that the first sample does
//              not start from t = 0.0 by incorporating a time offset.

#ifndef included_NonuniformDMD_h
#define included_NonuniformDMD_h

#include "DMD.h"

namespace CAROM {

/**
 * Class NonuniformDMD implements the NonuniformDMD algorithm on a given snapshot matrix.
 */
class NonuniformDMD : public DMD
{
public:

    /**
     * @brief Constructor.
     *
     * @param[in] dim        The full-order state dimension.
     */
    NonuniformDMD(int dim);

private:

    /**
     * @brief Unimplemented default constructor.
     */
    NonuniformDMD();

    /**
     * @brief Unimplemented copy constructor.
     */
    NonuniformDMD(
        const NonuniformDMD& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    NonuniformDMD&
    operator = (
        const NonuniformDMD& rhs);

    /**
     * @brief Construct f_snapshots_minus and f_snapshots_plus
     */
    std::pair<Matrix*, Matrix*> computePlusMinusSnapshotMatrices(const Matrix* snapshots);

    /**
     * @brief Compute phi.
     */
    void computePhi(struct DMDInternal dmd_internal_obj);

    /**
     * @brief Compute the appropriate exponential function when predicting the solution.
     */
    std::complex<double> computeEigExp(std::complex<double> eig, double t);

};

}

#endif
