/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

/*!
 * @file ScalaMat.hpp
 * 
 * @brief This header defines the data structure used to work with all ScaLAPACK
 * matrices with its public interface.
 * 
 * @invariant The `m_data` private member is always a valid
 * pointer as this streamlines logic when passing it to ScaLAPACK, which will
 * not appropriately handle a NULL pointer even when it shouldn't need to be
 * accessed.
 * 
 * @invariant `m_ctxt` either points to a valid BLACS context or it is nullptr.
 */

#ifndef SCALAMAT_HPP
#define SCALAMAT_HPP

#include "Context.hpp"

#include <array>
#include <iostream>
#include <memory>

namespace ScalaWRAP {

class LocalMatrix;

///
/// If no user-specified blocksize is given, ScaLAPACK matrices are distributed
/// in 32 x 32 blocks.
///
constexpr static int default_blocksize = 32;

/*!
 * @brief The class encapsulating a distributed matrix for all operations.
 * Most operations will require all processes to call them, or else they will
 * hang.
 */
class ScalaMat
{
private:
    int m_globm, m_globn, m_mb, m_nb, m_rowsrc, m_colsrc;
    std::shared_ptr<const Context> m_ctxt;
    double* m_data;
    bool m_own_data;

    enum class Direction { ROW, COL };

    template<Direction D>
    struct TypedRange
    {
        int lower; int upper;
        TypedRange(int i, int j) : lower(i), upper(j) {}
    };

public:
    /*!
     * @brief This type comes up a lot in code using this class; have a
     * convenience typedef.
     */
    typedef std::shared_ptr<const Context> context_ptr;

    context_ptr context() const noexcept { return m_ctxt; }

    /*!
     * @brief Construct a ScaLAPACK matrix with given size (optionally
     * specifying BLACS parameters), but uninitialized data.
     * 
     * This constructor is the one you should use most often - the default
     * choices are good for most use cases and will mostly be manipulated
     * internally. The only one you may wish to adjust is the blocksize, if
     * there is an advantage to a different storage layout in your application.
     * 
     * @param[in] global_m The global number of rows of the distributed matrix
     * 
     * @param[in] global_n The global number of columns of the distributed
     * matrix
     * 
     * @param[in] mb (Default: default_blocksize) The blocking factor in the row
     * direction.
     * 
     * @param[in] nb (Default: default_blocksize) The blocking factor in the
     * column direction.
     * 
     * @param[in] ctxt (Default: Context::default_context()) The BLACS context
     * in which the matrix is distributed.
     * 
     * @param[in] rowsrc (Default: 0) The starting index for distribution of the
     * matrix in the row direction.
     * 
     * @param[in] colsrc (Default: 0) The starting index for distribution of the
     * matrix in the column direction.
     * 
     * @note For some ScaLAPACK arithmetic routines, `mb == nb` is a
     * prerequisite. There is also probably little reason for a user to change
     * `rowsrc`, `colsrc`, or `ctxt`.
     */
    ScalaMat(int global_m, int global_n, int mb = default_blocksize,
             int nb = default_blocksize,
             std::shared_ptr<const Context> ctxt = Context::default_context(),
             int rowsrc = 0, int colsrc = 0);

    /*!
     * @brief Constructor for expert use only.
     * 
     * This constructor is for use if you *really* need to know exactly where
     * your data is allocated and deallocated, re-use other storage, etc. You
     * are responsible for making sure that `data` is allocated with the
     * appropriate amount of storage on each process; you can determine the
     * dimensions using `ctxt->getinfo()` and `numroc`. If you think you need
     * to use this constructor, you should probably make sure to read the
     * ScaLAPACK user's guide.
     * ([SLUG](http://www.netlib.org/scalapack/slug/node1.html)).
     * 
     * The constructed object does __not__ take ownership of the data array;
     * you will be responsible for deallocating it.
     * 
     * @param[inout] data Pointer to memory allocated with at least the amount
     * of storage required to hold a ScaLAPACK matrix distributed with the
     * given parameters. This data is *not* copied and the constructed object
     * does *not* take ownership - you are responsible for deallocating it
     * later.
     * 
     * @param[in] global_m Row dimension of the distributed array.
     * 
     * @param[in] global_n Column dimension of the distributed array.
     * 
     * @param[in] mb The row blocking factor to be used in the distributed
     * array.
     * 
     * @param[in] nb The column blocking factor to be used in the distributed
     * array.
     * 
     * @param[in] ctxt Shared pointer to the context where this matrix lives.
     * The constructor does a sanity check to make sure that `rowsrc` and
     * `colsrc` are valid indices in the given context's process grid.
     * 
     * @param[in] rowsrc The 0-based row index of the process in the process
     * grid that distribution starts from.
     * 
     * @param[in] colsrc The 0-based column index of the process in the process
     * grid that distribution starts from.
     */
    ScalaMat(double* data, int global_m, int global_n, int mb, int nb,
             std::shared_ptr<const Context> ctxt, int rowsrc, int colsrc);

    /*!
     * @brief Copy construction is deleted to avoid surprises.
     * 
     * I didn't want to have either:
     *   1) Implicit copying of a large parallel matrix, or
     *   2) Implicit sharing of memory so that modifications to one matrix
     *      are reflected in another.
     * 
     * Therefore the copy assignment operator is deleted.
     * 
     * In the case that you want to have two `ScalaMat`'s with the same
     * underlying data, use `std::move`. The state of a moved from `ScalaMat` is
     * valid, with only the exception that it no longer owns its data and will
     * not delete it when the destructor is called. Example: `A = std::move(B)`.
     */
    ScalaMat(const ScalaMat& other) = delete;

    /*!
     * @brief Copy assignment operator is deleted for the same reason as copy
     * constructor. See detailed description on the copy constructor.
     */
    ScalaMat& operator=(const ScalaMat& other) = delete;

    /*!
     * @brief Construct a `ScalaMat` by taking ownership of data from `other`.
     * 
     * The result of this operation is that the constructed matrix now owns the
     * data from `other` if `other` owned it (otherwise neither matrix owns it,
     * wherever it came from is responsible for deletion) and `other` has no
     * ownership. `other` still represents a valid matrix, until such time as
     * the backing memory is freed.
     */
    ScalaMat(ScalaMat&& other);

    /*!
     * @brief Construct a localized `ScalaMat` whose data is only on a single
     * process. The LocalMatrix needs to have been initialized on all processes
     * for communication purposes.
     * 
     * `other` must be stored in column major order; this constructor will throw
     * an exception if it is not.
     * 
     * This constructor takes ownership of the data from `other` if `other` had
     * it to begin with. Make sure that memory is externally managed and `other`
     * does not own it if you intend to re-use it.
     */
    ScalaMat(LocalMatrix&& other,
             std::shared_ptr<const Context> ctxt = Context::default_context());

    /*!
     * @brief Copy the data from this matrix into `other` in transposed order.
     * `other` must be allocated with the appropriate size.
     */
    ScalaMat& transpose(ScalaMat& other) const;

    /*!
     * @brief Take ownership of matrix data from `other`. See note on the move
     * constructor.
     */
    ScalaMat& operator=(ScalaMat&& other);

    int m() const noexcept { return m_globm; }
    int n() const noexcept { return m_globn; }
    int mb() const noexcept { return m_mb; }
    int nb() const noexcept { return m_nb; }
    int rowsrc() const noexcept { return m_rowsrc; }
    int colsrc() const noexcept { return m_colsrc; }

    /*!
     * @brief Represents a view of an underlying distributed matrix.
     * 
     * I disallowed copying data between ScaLAPACK matrices because this is
     * probably not usually a good idea and we shouldn't be able to do it by
     * accident. However, we do want to be able to move data around. The
     * `Submatrix` class provides a view that doesn't actually own any data, and
     * can be used for assignment operations. `A = B` could then be achieved by
     *
     *     A.submatrix(rowrange(1, A.m()), colrange(1, A.n())) =
     *          B.submatrix(rowrange(1, B.m()), colrange(1, B.n()));
     * 
     * assuming, of course, that `A` and `B` have compatible size. This notation
     * is intentionally verbose because copying an entire parallel matrix is
     * not usually the best thing to do.
     */
    class Submatrix
    {
    public:
        /*!
         * @brief Assign the contents of one Submatrix to another. This copies
         * data between the underlying distributed matrices.
         */
        Submatrix& operator=(const Submatrix& other);

        /*!
         * @brief Assign the entire contents of a `LocalMatrix` to a Submatrix
         * of a distributed matrix. Use to distribute data from individual
         * processors to the distributed grid.
         * 
         * Because ScaLAPACK matrices are stored in column major order, the
         * assignment will perform a transpose operation on `other` if it is
         * stored in row major order (by default a `LocalMatrix` is assumed to
         * be row major when constructed). This operation is performed in the
         * constructor overload `ScalaMat(const LocalMatrix&)`.
         */
        Submatrix& operator=(LocalMatrix&& other);

        Submatrix& operator+=(const Submatrix& other);
        Submatrix& operator+=(LocalMatrix&& other);

        int m() const noexcept { return m_m; }
        int n() const noexcept { return m_n; }
        int i() const noexcept { return m_i; }
        int j() const noexcept { return m_j; }

        const ScalaMat& parent() const noexcept { return m_parent; }
        ScalaMat& parent() noexcept { return m_parent; }

        ScalaMat similar() const;

        void initialize_descriptor(std::array<int, 9>& d) const noexcept
        {
            m_parent.initialize_descriptor(d);
            d[2] = m();
            d[3] = n();
        }

    private:
        friend class ScalaMat;
        Submatrix(ScalaMat& parent, std::pair<int, int> ibounds, 
              std::pair<int, int> jbounds);

        ScalaMat& m_parent;
        int m_i, m_j, m_m, m_n;
    };

    /*!
     * @brief For clarity, when taking a view of the matrix with `Submatrix`,
     * the Submatrix bounds must be specified with these wrapper types. Use
     * `RowRange` to indicate, well, a range of rows.
     */
    typedef TypedRange<Direction::ROW> RowRange;

    /*!
     * @brief For clarity, when taking a view of the matrix with `Submatrix`,
     * the Submatrix bounds must be specified with these wrapper types. Use
     * `ColRange` to indicate a range of columns.
     */
    typedef TypedRange<Direction::COL> ColRange;

    /*!
     * @brief Get a "view" of the submatrix
     * `A[ibounds.lower:ibounds.upper, jbounds.lower:jbounds.upper]`
     * 
     * Indexing is 1-based, reflecting the fact that ScaLAPACK is a Fortran
     * library and also the fact that 1-based indexing is nice for linear
     * algebra. The returned `Submatrix` may be assigned to from another
     * `Submatrix`, or a `LocalMatrix`.
     */
    Submatrix submatrix(RowRange ibounds, ColRange jbounds);
    const Submatrix submatrix(RowRange ibounds, ColRange jbounds) const;

    Submatrix operator()(RowRange ibounds, ColRange jbounds)
    {
        return submatrix(ibounds, jbounds);
    }

    const Submatrix operator()(RowRange ibounds, ColRange jbounds) const
    {
        return submatrix(ibounds, jbounds);
    }

    ScalaMat similar() const
    {
        return ScalaMat(m(), n(), m_mb, m_nb, context(), m_rowsrc, m_colsrc);
    }

    const double* data() const noexcept { return m_data; }
    double* data() noexcept { return m_data; }
    double* release_data() noexcept { m_own_data = false; return m_data; }
    bool own_data() const noexcept { return m_own_data; }

    ScalaMat& operator=(double x) noexcept;

    ~ScalaMat()
    {
        if (m_own_data && m_data != nullptr) {
            delete[] m_data;
        }
    }

    void initialize_descriptor(std::array<int, 9>&) const noexcept;

    void dump_debug_info(std::ostream& o, bool show_data = false) const;

    static void print_descriptor(std::ostream&, const std::array<int, 9>&);
};

inline ScalaMat::RowRange rowrange(int i, int j)
{
    return ScalaMat::RowRange(i, j);
}

inline ScalaMat::ColRange colrange(int i, int j)
{
    return ScalaMat::ColRange(i, j);
}

} /* namespace ScalaWRAP */

#endif /* SCALAMAT_HPP */
