/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: An abstract class defining the interface to the generic SVD
//              algorithm.

#include "SVD.h"

namespace CAROM {

SVD::SVD(
    Options options) :
    d_dim(options.dim),
    d_num_samples(0),
    d_max_num_samples(options.max_num_samples),
    d_debug_algorithm(options.debug_algorithm)
{
    CAROM_VERIFY(options.dim > 0);
    CAROM_VERIFY(options.max_num_samples > 0);
}

}
