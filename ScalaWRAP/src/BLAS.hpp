/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

#ifndef SCALAWRAP_BLAS_HPP
#define SCALAWRAP_BLAS_HPP

/*!
 * @file BLAS.hpp
 * 
 * @brief Declares the interface to BLAS functionality for ScalaWRAP objects.
 * 
 * This file is mainly used inside the library for testing - the user may find
 * the PBLAS functionality wrapped here very useful, however.
 */

#include "ScalaMat.hpp"
#include "LocalMatrix.hpp"

namespace ScalaWRAP {

using Submatrix = ScalaMat::Submatrix;

Submatrix& axpby(bool trans, double alpha, const Submatrix& X, double beta,
                 Submatrix& Y);
ScalaMat& axpby(bool trans, double alpha, const ScalaMat& X, double beta, ScalaMat& Y);
LocalMatrix& axpby(bool trans, double alpha, const LocalMatrix& X, double beta,
                   LocalMatrix& Y);

inline Submatrix& axpy(bool trans, double alpha, const Submatrix& X, Submatrix& Y)
{
    return axpby(trans, alpha, X, 1.0, Y);
}

inline ScalaMat& axpy(bool trans, double alpha, const ScalaMat& X, ScalaMat& Y)
{
    return axpby(trans, alpha, X, 1.0, Y);
}

inline LocalMatrix& axpy(bool trans, double alpha, const LocalMatrix& X, LocalMatrix& Y)
{
    return axpby(trans, alpha, X, 1.0, Y);
}

double norm(const LocalMatrix& x);

ScalaMat& gemm(bool transa, double alpha, const ScalaMat& A,
               bool transb, const ScalaMat& B, double beta, ScalaMat& C);

Submatrix& gemm(bool transa, double alpha, const Submatrix& A,
                bool transb, const Submatrix& B, double beta, Submatrix& C);

} /* namespace ScalaWRAP */

#endif /* SCALAWRAP_BLAS_HPP */