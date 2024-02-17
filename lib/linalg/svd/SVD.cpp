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
    d_max_num_samples(options.samples_per_time_interval),
    d_basis(NULL),
    d_basis_right(NULL),
    d_U(NULL),
    d_W(NULL),
    d_S(NULL),
    d_snapshots(NULL),
    d_debug_algorithm(options.debug_algorithm)
{
    CAROM_VERIFY(options.dim > 0);
    CAROM_VERIFY(options.samples_per_time_interval > 0);
}

SVD::~SVD()
{
    delete d_basis;
    delete d_U;
    delete d_S;
    delete d_basis_right;
    delete d_W;
    delete d_snapshots;
}

}
