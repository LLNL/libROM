/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

/*!
 * @file LocalMatrix.hpp
 * 
 * @brief This header defines the `LocalMatrix` class used as a go-between for
 * local data and distributed matrices encapsulated by `ScalaMAT`.
 */
#ifndef LOCALMATRIX_HPP
#define LOCALMATRIX_HPP

#include <exception>

#include "MPI_utils.hpp"
#include "ScalaMat.hpp"

namespace ScalaWRAP {

struct LocalMatrixException : public std::exception
{
    const char* what() const noexcept
    {
        return "m_data is nullptr but m != 0 and n != 0";
    }
};

/*!
 * @brief Sentinel values to indicate the memory ordering of a `LocalMatrix`.
 * By default, local arrays are assumed to be row major.
 */
enum class MemoryOrder
{
    ROW_MAJOR, COL_MAJOR
};

/// Export the ROW_MAJOR constant to be used by other routines
constexpr auto ROW_MAJOR = MemoryOrder::ROW_MAJOR;

/// Export the COL_MAJOR constant to be used by other routines
constexpr auto COL_MAJOR = MemoryOrder::COL_MAJOR;

/*!
 * @brief Class encapsulating a matrix (2-dimensional array) stored on a single
 * process.
 * 
 * ScalaWRAP does not make great use of local matrices since its purpose is to
 * wrap ScaLAPACK operations on distributed matrices. It does not re-export
 * BLAS or LAPACK functionality to deal with operations on localized data. There
 * is a need, however, for an object to act as a go-between between the local
 * and distributed scopes, which is what the `LocalMatrix` class provides.
 * 
 * A `LocalMatrix` contains a pointer to some data, which it may or may not own.
 * If owned, the object deletes the data when destroyed if it hasn't passed
 * ownership to another object. By default data is assumed not to be owned. The
 * matrix has a size `m x n` specified by the user, a member tracking which
 * process actually owns the data (the data pointer on other processes will not
 * be accessed), and a flag to indicate whether the data is in C or Fortran
 * order (row vs. column major).
 * 
 * The only meaningful operation provided on `LocalMatrix`'s is using them to
 * collect or broadcast data to a global `ScalaMat`. A `LocalMatrix` may be
 * assigned to any submatrix of a `ScalaMat`, and any `Submatrix` may be
 * collected to a `LocalMatrix` of appropriate size using `operator=`. As a
 * special case, a `ScalaMat` may be constructed from a _column major_
 * `LocalMatrix`; the new object is a distributed matrix whose data is entirely
 * stored on a single process.
 */
class LocalMatrix
{
public:
    /*!
     * @brief Data constructor. Given parameters of a local array, construct a
     * `LocalMatrix`.
     * 
     * @param[inout] data Pointer to local array. If `take_data == true`, the
     * new object takes control of this array and may delete it or pass
     * ownership to another object.
     * 
     * @param[in] m The number of rows of the array.
     * 
     * @param[in] n The number of columns of the array.
     * 
     * @param[in] srcrank The rank of the process (in MPI_COMM_WORLD) on which
     * the local data resides. `data` is not referenced by other processes and
     * can point wherever.
     * 
     * @param[in] take_data (Default: false) Whether or not the new object owns
     * the data array. By default the data is owned by the caller and will not
     * be freed - if ownership is given to the `LocalMatrix` it will be
     * difficult to track where the data is actually freed and whether it is
     * safe to access. This is because, for performance reasons, the ownership
     * may be handed off to a `ScalaMat` object in a subsequent call, and that
     * object will then deallocate the data itself. As a recommendation,
     * consider allocating local arrays using a smart pointer and allowing the
     * pointer to free the data when you are done with it, and leave `take_data`
     * as `false`.
     * 
     * @param[in] order The ordering of the array; by default, row major
     * ordering is assumed. Scatter and gather operations take the ordering into
     * account.
     */
    LocalMatrix(double* data, int m, int n, int srcrank, bool take_data = false,
                MemoryOrder order = ROW_MAJOR)
        : m_data(data), m_m(m), m_n(n), m_rank(srcrank), m_ordering(order), 
          m_own_data(take_data)
    {
        if (mpi_rank() != srcrank) {
            if (m_data && m_own_data) { delete[] m_data; }
            m_data = nullptr; m_own_data = false;
        } else {
            if (m_data == nullptr && m != 0 && n != 0) {
                throw LocalMatrixException();
            }
        }
    }

    /*!
     * @brief Default constructor makes a local matrix that is not associated
     * to any process and has no data.
     * 
     * @todo Is this actually used somewhere?
     */
    LocalMatrix() : m_data(nullptr), m_m(0), m_n(0), m_rank(-1),
                    m_ordering(ROW_MAJOR), m_own_data(false)
    {}

    ~LocalMatrix()
    {
        if (m_own_data && m_data != nullptr) {
            delete[] m_data;
        }
    }

    /*!
     * @brief Copy construction and assignment is disabled.
     * 
     * Copy semantics don't make sense with the architecture of the class.
     * Move construction leaves the source object in a valid state, just without
     * ownership of its data array.
     */
    LocalMatrix(const LocalMatrix& other) = delete;
    LocalMatrix& operator=(const LocalMatrix& other) = delete;

    /*!
     * @brief Move construction takes ownership of the data from `other` and
     * duplicates all other parameters.
     * 
     * `other` is left as a valid `LocalMatrix`, but without ownership of its
     * data it may be difficult to know when the pointer is freed. If `other`
     * did not own its data array, the newly constructed object does not either
     * and it is still safe to use.
     */
    LocalMatrix(LocalMatrix&& other) : m_data(other.m_data), m_m(other.m_m),
                                       m_n(other.m_n), m_rank(other.m_rank),
                                       m_ordering(other.m_ordering),
                                       m_own_data(other.m_own_data)
    {
        other.m_own_data = false;
    }

    /*!
     * @brief Move assignment is supported with the same semantics as the move
     * constructor.
     */
    LocalMatrix& operator=(LocalMatrix&& other)
    {
        if (m_own_data && m_data != nullptr) {
            delete[] m_data;
        }
        new(this) LocalMatrix(std::move(other));
        return *this;
    }

    /*!
     * @brief Copy the contents of a submatrix of a distributed matrix into the
     * local data array.
     * 
     * This operation does account for memory ordering. In one special case, an
     * extra matrix must be allocated: if the local matrix is on a process that
     * is not a member of `other.parent().context()`, the global matrix must be
     * copied to a matrix on the local process, then that matrix is transposed
     * to reverse the memory ordering.
     * 
     * This operation should throw an exception if sizes are incompatible, but
     * there is at least one code path that will make it to ScaLAPACK code and
     * result in an error message from PXERBLA and an abort instead of catching
     * the incompatible sizes.
     * 
     * @todo Make sure all error cases are covered by exceptions.
     */
    LocalMatrix& operator=(const ScalaMat::Submatrix& other);

    int srcrank() const noexcept { return m_rank; }
    MemoryOrder ordering() const noexcept { return m_ordering; }
    int m() const noexcept { return m_m; }
    int n() const noexcept { return m_n; }

    /// Return a const pointer to my data for observing.
    const double* data() const noexcept { return m_data; }

    /// Return a mutable pointer to my data, for passing to external routines.
    double* release_data() noexcept
    {
        m_own_data = false;
        return m_data;
    }

    bool own_data() const noexcept { return m_own_data; }
    
private:
    double* m_data;
    int m_m, m_n, m_rank;
    MemoryOrder m_ordering;
    bool m_own_data;
};

} /* namespace ScalaWRAP */

#endif /* LOCALMATRIX_HPP */
