/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

#include "BLAS.hpp"
#include "C_interface.h"

#include <exception>

namespace ScalaWRAP
{

struct BLASException : public std::exception
{
    const char *msg;
    const char *what() const noexcept { return msg; }
    BLASException(const char *m) : msg(m) {}
};

Submatrix &axpby(bool trans, double alpha, const Submatrix &X, double beta,
                 Submatrix &Y)
{
    if (!X.parent().context())
    {
        return Y;
    }
    if (X.parent().context()->id() != Y.parent().context()->id())
    {
        throw BLASException("Context mismatch in axpby (matrices must share a context");
    }
    else if (trans)
    {
        if (X.m() != Y.n() || X.n() != Y.m())
        {
            throw BLASException("Dimension mismatch in axpby");
        }
    }
    else
    {
        if (X.m() != Y.m() || X.n() != Y.n())
        {
            throw BLASException("Dimension mismatch in axpby");
        }
    }

    std::array<int, 9> descx, descy;
    X.initialize_descriptor(descx);
    Y.initialize_descriptor(descy);

    pdgeadd_wrapper(trans ? 'T' : 'N', Y.m(), Y.n(), alpha, X.parent().data(),
                    X.i(), X.j(), descx.data(), beta, Y.parent().data(), Y.i(),
                    Y.j(), descy.data());
    return Y;
}

ScalaMat &axpby(bool trans, double alpha, const ScalaMat &X, double beta,
                ScalaMat &Y)
{
    auto Xsub = X(rowrange(1, X.m()), colrange(1, X.n()));
    auto Ysub = Y(rowrange(1, Y.m()), colrange(1, Y.n()));
    axpby(trans, alpha, Xsub, beta, Ysub);
    return Y;
}

LocalMatrix &axpby(bool trans, double alpha, const LocalMatrix &X, double beta,
                   LocalMatrix &Y)
{
    if (trans)
    {
        if (X.m() != Y.n() || X.n() != Y.m())
        {
            throw BLASException("Dimension mismatch in axpby");
        }
    }
    else
    {
        if (X.m() != Y.m() || X.n() != Y.n())
        {
            throw BLASException("Dimension mismatch in axpby");
        }
    }

    bool owned = Y.own_data();
    double *Ydata = Y.release_data();
    axpby_wrapper(X.m() * X.n(), alpha, X.data(), 1, beta, Ydata, 1);
    Y = LocalMatrix(Ydata, Y.m(), Y.n(), Y.srcrank(), owned, Y.ordering());
    return Y;
}

double norm(const LocalMatrix &x)
{
    return dnrm2_wrapper(x.m() * x.n(), x.data(), 1);
}

Submatrix &gemm(bool transa, double alpha, const Submatrix &A, bool transb,
                const Submatrix &B, double beta, Submatrix &C)
{
    if (!A.parent().context())
    {
        return C;
    }

    // Check that contexts match
    if (A.parent().context() != B.parent().context() ||
        B.parent().context() != C.parent().context())
    {
        throw BLASException("Matrix contexts must match");
    }

    // Check dimensions.
    if (transa) {
        if (transb) {
            if (A.m() != B.n() || A.n() != C.m() || B.m() != C.n()) {
                throw BLASException("Dimension mismatch in gemm");
            }
        } else {
            if (A.m() != B.m() || A.n() != C.m() || B.n() != C.n()) {
                throw BLASException("Dimension mismatch in gemm");
            }
        }
    } else {
        if (transb) {
            if (A.n() != B.n() || A.m() != C.m() || B.m() != C.n()) {
                throw BLASException("Dimension mismatch in gemm");
            }
        } else {
            if (A.n() != B.m() || A.m() != C.m() || B.n() != C.n()) {
                throw BLASException("Dimension mismatch in gemm");
            }
        }
    }

    std::array<int, 9> desca, descb, descc;
    A.initialize_descriptor(desca);
    B.initialize_descriptor(descb);
    C.initialize_descriptor(descc);

    pdgemm_wrapper(transa ? 'T' : 'N', transb ? 'T' : 'N', C.m(), C.n(),
                   transa ? A.m() : A.n(), alpha, A.parent().data(), A.i(),
                   A.j(), desca.data(), B.parent().data(), B.i(), B.j(),
                   descb.data(), beta, C.parent().data(), C.i(), C.j(),
                   descc.data());
    return C;
}

ScalaMat& gemm(bool transa, double alpha, const ScalaMat& A, bool transb,
               const ScalaMat& B, double beta, ScalaMat& C)
{
    auto Asub = A(rowrange(1, A.m()), colrange(1, A.n()));
    auto Bsub = B(rowrange(1, B.m()), colrange(1, B.n()));
    auto Csub = C(rowrange(1, C.m()), colrange(1, C.n()));

    gemm(transa, alpha, Asub, transb, Bsub, beta, Csub);
    return C;
}

} /* namespace ScalaWRAP */