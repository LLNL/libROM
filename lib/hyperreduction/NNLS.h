/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Interface to the NNLS algorithm for sparse approximation of a
// non-negative least squares (NNLS) solution.

#ifndef included_NNLS_h
#define included_NNLS_h

#include "linalg/Matrix.h"
#include "linalg/Vector.h"

namespace CAROM {

/**
 * @brief Computes the approximate NNLS solution x for a given system matrix G,
 * tolerance tau, and RHS Gw.
 *
 * Implemented from T. Chapman et al, "Accelerated mesh sampling for the hyper
 * reduction of nonlinear computational models," Int. J. Numer. Meth. Engng.,
 * 109: 1623-1654.
 *
 * @param[in] tau The positive tolerance used for termination of iterations.
 * @param[in] G The dense system matrix.
 * @param[in] w The exact dense solution, which defines the RHS as Gw.
 * @param[out] x The sparse approximate NNLS solution.
 * @param[in] printLevel If nonzero, iteration output is written.
 */
void NNLS(const double tau,
	  Matrix const& G,
	  Vector const& w,
	  Vector & x,
	  int printLevel = 0);

}

#endif
