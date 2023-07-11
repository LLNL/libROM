/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: A class implementing interface of SVD for the Randomized SVD
//              algorithm.

#ifndef included_RandomizedSVD_h
#define included_RandomizedSVD_h

#include "StaticSVD.h"
#include "linalg/Options.h"
#include "linalg/scalapack_wrapper.h"

#include <limits>
#include <memory>
#include <vector>

namespace CAROM {

/**
 * Class RandomizedSVD implements the Randomized SVD algorithm introduced in
 *    Nathan Halko, Per-Gunnar Martinsson, and Joel A. Tropp.
 *    "Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions."
 *    SIAM review 53.2 (2011): 217-288.
 */
class RandomizedSVD : public StaticSVD
{
private:
    friend class BasisGenerator;

    /**
     * @brief Constructor.
     *
     * @param[in] options The struct containing the options for this SVD
     *                    implementation.
     */
    RandomizedSVD(
        Options options
    );

    /**
     * @brief Unimplemented default constructor.
     */
    RandomizedSVD();

    /**
     * @brief Unimplemented copy constructor.
     */
    RandomizedSVD(
        const RandomizedSVD& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    RandomizedSVD&
    operator = (
        const RandomizedSVD& rhs);

    /**
     * @brief Gathers samples from all other processors to form complete
     * sample of system and computes the SVD.
     */
    void
    computeSVD();

    /**
     * @brief The number of dimensions of the randomized subspace the
     * snapshot matrix will be projected to.
     */
    int d_subspace_dim;
};

}

#endif
