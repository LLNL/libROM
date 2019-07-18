/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

#include "SVD.hpp"
#include "C_interface.h"
#include "MPI_utils.hpp"

#include <algorithm>
#include <iostream>

namespace ScalaWRAP
{

SVDInfo svd(ScalaMat::Submatrix &A, bool wantu, bool wantv)
{
    if (!A.parent().context())
    {
        return SVDInfo();
    }

    // Follow naming from the ScaLAPACK doc comments.
    const int SIZE = std::min(A.m(), A.n());

    SVDInfo info = {nullptr, nullptr, {}, 0};
    info.S.resize(SIZE);

    if (wantu)
    {
        info.U = std::unique_ptr<ScalaMat>(
            new ScalaMat(A.m(), SIZE, A.parent().mb(), A.parent().nb(),
                         A.parent().context(), A.parent().rowsrc(),
                         A.parent().colsrc()));
    }

    if (wantv)
    {
        info.Vt = std::unique_ptr<ScalaMat>(
            new ScalaMat(SIZE, A.n(), A.parent().mb(), A.parent().nb(),
                         A.parent().context(), A.parent().rowsrc(),
                         A.parent().colsrc()));
    }

    const char jobu = wantu ? 'V' : 'N';
    const char jobvt = wantv ? 'V' : 'N';
    static double dummy = 0;
    auto Adata = A.parent().data();
    auto Udata = info.U ? info.U->data() : &dummy;
    auto Vdata = info.Vt ? info.Vt->data() : &dummy;
    std::array<int, 9> desca, descu, descvt;
    A.initialize_descriptor(desca);
    if (info.U)
    {
        info.U->initialize_descriptor(descu);
    }
    if (info.Vt)
    {
        info.Vt->initialize_descriptor(descvt);
    }

    std::vector<double> work{0.0};

    info.info = pdgesvd_wrapper(
        jobu, jobvt, A.m(), A.n(), Adata, A.i(), A.j(), desca.data(),
        info.S.data(), Udata, 1, 1, descu.data(), Vdata, 1, 1, descvt.data(),
        work.data(), -1);

    work.resize(static_cast<unsigned>(work[0]));

    info.info = pdgesvd_wrapper(
        jobu, jobvt, A.m(), A.n(), Adata, A.i(), A.j(), desca.data(),
        info.S.data(), Udata, 1, 1, descu.data(), Vdata, 1, 1, descvt.data(),
        work.data(), work.size());

    if (info.info < 0)
    {
        if (mpi_rank() == 0)
        {
            std::cerr << "Error in PDGESVD - invalid argument; return value was "
                      << info.info << '\n';
            std::cerr << "(in " << __func__ << " at " << __FILE__ << ':'
                      << __LINE__ << ")\n";
        }
        throw SVDException("Invalid argument");
    }
    else if (info.info > 0)
    {
        if (mpi_rank() == 0)
        {
            if (info.info == SIZE + 1)
            {
                std::cerr << "Warning; singular values not indentical across processes\n";
            }
            else
            {
                std::cerr << "Warning: DBDSQR did not converge in SVD\n";
            }
            std::cerr << "(in " << __func__ << " at " << __FILE__ << ':'
                      << __LINE__ << ")\n";
        }
    }
    return info;
}

} /* namespace ScalaWRAP */