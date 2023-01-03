/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

#include "Utilities.hpp"

namespace CAROM {

void ComputeCtAB(const HypreParMatrix& A,
                 const CAROM::Matrix& B,  // Distributed matrix
                 const CAROM::Matrix& C,  // Distributed matrix
                 CAROM::Matrix& CtAB)     // Non-distributed (local) matrix
{
    MFEM_VERIFY(B.distributed() && C.distributed() && !CtAB.distributed(),
                "In ComputeCtAB, B and C must be distributed, but not CtAB.");

    const int num_rows = B.numRows();
    const int num_cols = B.numColumns();
    const int num_rows_A = A.GetNumRows();

    MFEM_VERIFY(C.numRows() == num_rows_A, "");

    mfem::Vector Bvec(num_rows);
    mfem::Vector ABvec(num_rows_A);

    CAROM::Matrix AB(num_rows_A, num_cols, true);

    for (int i = 0; i < num_cols; ++i) {
        for (int j = 0; j < num_rows; ++j) {
            Bvec[j] = B(j, i);
        }
        A.Mult(Bvec, ABvec);
        for (int j = 0; j < num_rows_A; ++j) {
            AB(j, i) = ABvec[j];
        }
    }

    C.transposeMult(AB, CtAB);
}

void ComputeCtAB_vec(const HypreParMatrix& A,
                     const HypreParVector& B, // Distributed vector
                     const CAROM::Matrix& C,  // Distributed matrix
                     CAROM::Vector& CtAB_vec) // Non-distributed (local) vector
{
    MFEM_VERIFY(C.distributed() && !CtAB_vec.distributed(),
                "In ComputeCtAB_vec, C must be distributed, but not CtAB_vec");

    MFEM_VERIFY(C.numRows() == A.GetNumRows(), "");
    MFEM_VERIFY(B.GlobalSize() == A.GetGlobalNumRows(), "");

    HypreParVector* AB = new HypreParVector(B);
    A.Mult(B, *AB);

    CAROM::Vector AB_carom(AB->GetData(), AB->Size(), true);
    C.transposeMult(AB_carom, CtAB_vec);
}

} // end namespace CAROM
