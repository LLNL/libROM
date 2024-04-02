/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
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
    const int num_rows_A = A.NumRows();

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

    MFEM_VERIFY(C.numRows() == A.NumRows(), "");
    MFEM_VERIFY(B.GlobalSize() == A.GetGlobalNumRows(), "");

    HypreParVector* AB = new HypreParVector(B);
    A.Mult(B, *AB);

    CAROM::Vector AB_carom(AB->GetData(), AB->Size(), true);
    C.transposeMult(AB_carom, CtAB_vec);
}

void verify_within_portion(const mfem::Vector &bb_min,
                           const mfem::Vector &bb_max,
                           const mfem::Vector &t, const double limit)
{
    // helper function to check if t is within limit percentage relative
    //  to the center of the mesh
    CAROM_VERIFY(t.Size() == bb_min.Size() && bb_min.Size() == bb_max.Size());
    CAROM_VERIFY(limit >= 0.0 && limit <= 1.0);
    for (int i = 0; i < t.Size(); i++)
    {
        double domain_limit = limit * (bb_max[i] - bb_min[i]);
        double mesh_center = 0.5 * (bb_max[i] + bb_min[i]);

        // check that t is within the limit relative to the center of the mesh
        if (std::abs((t(i) - mesh_center)) - (0.5 * domain_limit) > 1.0e-14)
        {
            std::cerr << "Error: value of t exceeds domain limit: t = " << t(
                          i) << ", limit = " << 0.5 * domain_limit << "\n";
            exit(-1);
        }
    }
}

double map_to_ref_mesh(const double &bb_min, const double &bb_max,
                       const double &fraction)
{
    // helper function to map a fractional value from [-1, 1] to [bb_min, bb_max]
    CAROM_VERIFY(fraction <= 1.0 && fraction >= -1.0);
    return bb_min + (fraction + 1.0) * ((bb_max - bb_min) * 0.5);
}

double map_from_ref_mesh(const double &bb_min, const double &bb_max,
                         const double &value)
{
    // helper function to map a value from the mesh range [bb_min, bb_max] to [-1, 1]
    CAROM_VERIFY(value <= bb_max && value >= bb_min);
    return -1.0 + (value - bb_min) * ((2.0) / (bb_max - bb_min));
}

} // end namespace CAROM
