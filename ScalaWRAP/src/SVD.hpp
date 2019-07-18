/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

/*!
 * @file SVD.hpp
 * 
 * @brief This header declares the interface to a distributed SVD using
 * ScaLAPACK.
 * 
 * The only function exported here is `svd`, along with the exception type
 * that it throws if there is an argument error, and the info struct that it
 * returns.
 */

#ifndef SCALAWRAP_SVD_HPP
#define SCALAWRAP_SVD_HPP

#include "ScalaMat.hpp"

#include <exception>

namespace ScalaWRAP {

/*!
 * @brief Exception type thrown by svd
 * 
 * This exception is only thrown in case of an invalid argument to PDGESVD,
 * which should never happen.
 */
struct SVDException : public std::exception
{
    const char* msg;
    const char* what() const noexcept { return msg; }
    SVDException(const char* c) : msg(c) {}
};

/*!
 * @brief Encapsulates the returned information from the SVD.
 */
struct SVDInfo
{
    /// Pointer to the computed U, or nullptr if `wantu` is specified as `false`
    std::unique_ptr<ScalaMat> U;

    /// Pointer to the computed `V^T`, or nullptr if `wantv` is `false`.
    std::unique_ptr<ScalaMat> Vt;

    /// The vector of singular values, sorted in descending order of magnitude
    std::vector<double> S;

    /// The return code from PDGESVD; if non-zero you should have seen an error message.
    int info;
};

/*!
 * @brief Compute the SVD of `A`.
 * 
 * Computes the singular value decomposition of a (submatrix of) a distributed
 * `A`. The results are returned in the `SVDInfo` struct. This routine should
 * only ever fail due to numerical issues, since every matrix has an SVD.
 * 
 * @param[inout] A A `Submatrix`. The SVD of `A` is computed and the data stored
 * in `A` is destroyed in the process.
 * 
 * @param[in] wantu (Default: `true`) Set to `false` if you don't need the
 * matrix `U` of left singular vectors.
 * 
 * @param[in] wantvt (Default: `true`) Set to `false` if you don't need the
 * matrix `Vt` of transposed right singular vectors.
 * 
 * @returns An object of type `SVDInfo`. This contains pointers to `U` and `Vt`,
 * if they were requested, and the singular values in a `std::vector` as well
 * as the return code from PDGESVD. If `A` is `m x n` and `SIZE = min(m, n)`,
 * `U` is `m x SIZE`, there are `SIZE` singular values, and `Vt` is
 * `SIZE x n`.
 */
SVDInfo svd(ScalaMat::Submatrix& A, bool wantu = true, bool wantvt = true);

/*!
 * @brief Compute the SVD of `A`.
 * 
 * This overload is a wrapper that forwards the computation to the overload
 * for submatrices with a submatrix that contains the whole of `A`.
 */
inline SVDInfo svd(ScalaMat& A, bool wantu = true, bool wantvt = true)
{
    auto Asub = A(rowrange(1, A.m()), colrange(1, A.n()));
    return svd(Asub, wantu, wantvt);
}

} /* namespace ScalaWRAP */

#endif /* SCALAWRAP_SVD_HPP */
