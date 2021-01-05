/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
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
#include "Options.h"
#include "scalapack_wrapper.h"

#include <limits>
#include <memory>
#include <vector>

namespace CAROM {

/**
 * RandomizedSVD implements the interface of class SVD for the Randomized SVD
 * algorithm.  This algorithm is not scalable and is intended primarily as a
 * sanity check of the incremental svd algorithm.
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
       * @brief Number of dimensions of the randomized subspace the
       * snapshot matrix will be projected to.
       */
      int d_subspace_dim;
};

}

#endif
