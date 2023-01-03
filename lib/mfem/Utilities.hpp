/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

#ifndef MFEMUTILITIES_H
#define MFEMUTILITIES_H

#include "mfem.hpp"
#include "linalg/Matrix.h"

using namespace mfem;
using namespace std;

namespace CAROM {

/**
 * @brief This function computes a reduced, non-distributed matrix C^t AB stored
 *        identically (redundantly) on every process.
 *
 * @param[in] A The distributed HypreParMatrix (an MFEM class) A in C^t AB.
 *
 * @param[in] B The distributed Matrix B in C^t AB.
 *
 * @param[in] C The distributed Matrix C in C^t AB.
 *
 * @param[out] CtAB The non-distributed Matrix C^t AB.
 *
 */
void ComputeCtAB(const HypreParMatrix& A,
                 const CAROM::Matrix& B,
                 const CAROM::Matrix& C,
                 CAROM::Matrix& CtAB);

/**
 * @brief This function computes a reduced, non-distributed vector C^t AB stored
 *        identically (redundantly) on every process.
 *
 * @param[in] A The distributed HypreParMatrix (an MFEM class) A in C^t AB.
 *
 * @param[in] B The distributed HypreParVector B in C^t AB.
 *
 * @param[in] C The distributed Matrix C in C^t AB.
 *
 * @param[out] CtAB The non-distributed Vector C^t AB.
 *
 */
void ComputeCtAB_vec(const HypreParMatrix& A,
                     const HypreParVector& B,
                     const CAROM::Matrix& C,
                     CAROM::Vector& CtAB_vec);

}  // namespace CAROM

#endif // MFEMUTILITIES_H
