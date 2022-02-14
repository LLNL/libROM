/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
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
 *    "Finding Structure with Randomness: Probabilistic Algorithms for
 *    Constructing Approximate Matrix Decompositions" by N. Halko, P. G.
 *    Martinsson, and J. A. Tropp
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
