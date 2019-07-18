/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

#include "LocalMatrix.hpp"
#include "ScalaMat.hpp"
#include "C_interface.h"
#include "MPI_utils.hpp"

#include <algorithm>
#include <cassert>
#include <exception>
#include <iostream>
#include <utility>

namespace ScalaWRAP
{

struct ScalaMatException : public std::exception
{
    const char *msg;
    const char *what() const noexcept
    {
        return msg;
    }

    ScalaMatException(const char *m) noexcept : msg(m) {}
};

ScalaMat::ScalaMat(int global_m, int global_n, int mb, int nb,
                   std::shared_ptr<const Context> ctxt,
                   int rowsrc, int colsrc) : m_globm(global_m),
                                             m_globn(global_n),
                                             m_mb(mb), m_nb(nb),
                                             m_rowsrc(rowsrc), m_colsrc(colsrc),
                                             m_ctxt(ctxt), m_own_data(true)
{
    int nprow, npcol, pi, pj, mm, mn;
    if (ctxt)
    {
        blacs_gridinfo_wrapper(ctxt->id(), &nprow, &npcol, &pi, &pj);
    }
    if (ctxt == nullptr || nprow == -1)
    { // We are not in the context
        mm = mn = 0;
        m_ctxt = nullptr;
        m_globm = m_globn = m_mb = m_nb = m_rowsrc = m_colsrc = 0;
    }
    else
    {
        mm = numroc_wrapper(global_m, mb, pi, rowsrc, nprow);
        mn = numroc_wrapper(global_n, nb, pj, colsrc, npcol);
    }
    m_data = new double[mm * mn];
}

ScalaMat::ScalaMat(double* data, int global_m, int global_n, int mb, int nb,
                   std::shared_ptr<const Context> ctxt, int rowsrc,
                   int colsrc) : m_globm(global_m), m_globn(global_n), m_mb(mb),
                                 m_nb(nb), m_rowsrc(rowsrc), m_colsrc(colsrc),
                                 m_ctxt(ctxt), m_data(data), m_own_data(false)
{
    int nprow, npcol, pi, pj;
    std::tie(nprow, npcol, pi, pj) = ctxt->getinfo();
    if (rowsrc >= nprow || colsrc >= npcol) {
        std::cerr << "Error in expert constructor - given (rowsrc, colsrc) is "
                  << '(' << rowsrc << ',' << colsrc << ")\n, but the process grid "
                  << "is only " << nprow << " x " << npcol << '\n';
        std::cerr << "(in " << __func__ << " at " << __FILE__ << ':' << __LINE__ << ")\n";
        throw ScalaMatException("Specified rowsrc and colsrc outside process grid");
    }
}

ScalaMat &ScalaMat::operator=(double x) noexcept
{
    int nprow, npcol, pi, pj, mm, mn;
    if (m_ctxt)
    {
        std::tie(nprow, npcol, pi, pj) = m_ctxt->getinfo();
        mm = numroc_wrapper(m_globm, m_mb, pi, m_rowsrc, nprow);
        mn = numroc_wrapper(m_globn, m_nb, pj, m_colsrc, npcol);
        std::fill(m_data, m_data + size_t(mm) * mn, x);
    }
    return *this;
}

ScalaMat &ScalaMat::transpose(ScalaMat &other) const
{
    std::array<int, 9> my_desc, other_desc;
    initialize_descriptor(my_desc);
    other.initialize_descriptor(other_desc);

    pdgeadd_wrapper('T', other.m(), other.n(), 1.0, m_data, 1, 1, my_desc.data(),
                    0.0, other.m_data, 1, 1, other_desc.data());
    return other;
}

ScalaMat::ScalaMat(LocalMatrix &&other,
                   std::shared_ptr<const Context> ctxt) : m_globm(other.m()), m_globn(other.n()),
                                                          m_mb(std::max(other.m(), other.n())),
                                                          m_nb(m_mb), m_ctxt(ctxt), m_own_data(other.own_data())
{
    // Require that other is column major.
    if (other.ordering() == ROW_MAJOR)
    {
        if (mpi_rank() == 0)
        {
            std::cerr << "ScalaMat(LocalMatrix): Attempted to create a ScalaMat from a row major local\n";
            std::cerr << "(in " << __func__ << " at " << __FILE__ << ':'
                      << __LINE__ << ")\n";
        }
        throw ScalaMatException("It is not allowed to construct a ScalaMat from a row major matrix");
    }

    int nprow, npcol, pi = 0, pj = 0;

    if (ctxt)
    {
        blacs_gridinfo_wrapper(ctxt->id(), &nprow, &npcol, &pi, &pj);
    }

    int coords[] = {pi, pj};

    if (ctxt == nullptr || nprow == -1)
    { // We are not in the context.
        if (mpi_rank() == other.srcrank()) {
            std::cerr << "ScalaMat(LocalMatrix): Attempt to construct from a local matrix\n"
                      << "on a process outside the context\n";
            throw ScalaMatException("Local matrix not in context");
        }
        m_mb = m_nb = 0;
        m_rowsrc = m_colsrc = 0;
        m_ctxt = nullptr;
        m_data = new double[0];
        m_own_data = true;
        MPI_Bcast(coords, 2, MPI_INT, other.srcrank(), MPI_COMM_WORLD);
        return;
    }

    // Only processes in the participating context are still executing.
    // Set the coordinates of the source process in the grid.
    MPI_Bcast(coords, 2, MPI_INT, other.srcrank(), MPI_COMM_WORLD);
    m_rowsrc = coords[0];
    m_colsrc = coords[1];

    if (mpi_rank() == other.srcrank())
    {
        m_data = other.release_data();
    }
    else
    {
        m_data = new double[0];
        m_own_data = true;
    }
}

ScalaMat::ScalaMat(ScalaMat &&other) : m_globm(other.m_globm),
                                       m_globn(other.m_globn),
                                       m_mb(other.m_mb), m_nb(other.m_nb),
                                       m_rowsrc(other.m_rowsrc),
                                       m_colsrc(other.m_colsrc),
                                       m_ctxt(other.m_ctxt),
                                       m_data(other.m_data),
                                       m_own_data(other.m_own_data)
{
    other.m_own_data = false;
}

ScalaMat &ScalaMat::operator=(ScalaMat &&other)
{
    if (m_own_data && m_data != nullptr)
    {
        delete[] m_data;
    }
    new (this) ScalaMat(std::move(other));
    return *this;
}

ScalaMat::Submatrix &ScalaMat::Submatrix::operator=(
    const ScalaMat::Submatrix &other)
{
    int rank = mpi_global_params().first;
    if (m() != other.m() || n() != other.n())
    {
        if (rank == 0)
        {
            std::cerr << "Incompatible dimensions when assigning block\n"
                      << "dest is " << m() << " x " << n() << " but src is "
                      << other.m() << " x " << other.n() << '\n';
            std::cerr << "(in " << __func__ << " at " << __FILE__ << ':'
                      << __LINE__ << "\n";
        }
        throw ScalaMatException("Submatrix assignment: Incompatible dimensions");
    }

    std::array<int, 9> dst_descriptor, src_descriptor;
    m_parent.initialize_descriptor(dst_descriptor);
    other.m_parent.initialize_descriptor(src_descriptor);

    pdgemr2d_wrapper(m(), n(), other.m_parent.m_data, other.i(), other.j(),
                     src_descriptor.data(), m_parent.m_data, i(), j(),
                     dst_descriptor.data(), Context::default_context()->id());
    return *this;
}

ScalaMat::Submatrix &ScalaMat::Submatrix::operator+=(
    const ScalaMat::Submatrix &other)
{
    int rank = mpi_rank();
    if (m() != other.m() || n() != other.n())
    {
        if (rank == 0)
        {
            std::cerr << "Incompatible dimensions when assigning block\n"
                      << "dest is " << m() << " x " << n() << " but src is "
                      << other.m() << " x " << other.n() << '\n';
            std::cerr << "(in " << __func__ << " at " << __FILE__ << ':'
                      << __LINE__ << "\n";
        }
        throw ScalaMatException("Submatrix assignment: Incompatible dimensions");
    }
    else if (m_parent.context()->id() != other.parent().context()->id())
    {
        if (rank == 0)
        {
            std::cerr << "Matrices in arithmetic operations must share a context\n";
            std::cerr << "(in " << __func__ << " at " << __FILE__ << ':'
                      << __LINE__ << ")\n";
        }
        throw ScalaMatException("Submatrix addition: Mismatched contexts");
    }

    std::array<int, 9> dst_descriptor, src_descriptor;
    m_parent.initialize_descriptor(dst_descriptor);
    other.m_parent.initialize_descriptor(src_descriptor);

    pdgeadd_wrapper('N', m(), n(), 1.0, other.parent().data(), other.i(),
                    other.j(), src_descriptor.data(), 1.0, m_parent.m_data, i(),
                    j(), dst_descriptor.data());
    return *this;
}

ScalaMat::Submatrix &ScalaMat::Submatrix::operator+=(LocalMatrix &&other)
{
    int rank = mpi_rank();
    if (m() != other.m() || n() != other.n())
    {
        if (rank == 0)
        {
            std::cerr << "Incompatible dimensions when assigning block\n"
                      << "dest is " << m() << " x " << n() << " but src is "
                      << other.m() << " x " << other.n() << '\n';
            std::cerr << "(in " << __func__ << " at " << __FILE__ << ':'
                      << __LINE__ << "\n";
        }
        throw ScalaMatException("Submatrix assignment: Incompatible dimensions");
    }

    if (other.ordering() == COL_MAJOR)
    {
        ScalaMat A(std::move(other), parent().context());
        return (*this += A(rowrange(1, A.m()), colrange(1, A.n())));
    }
    else
    {
        bool owned = other.own_data();
        ScalaMat tmp1(LocalMatrix(other.release_data(), other.n(), other.m(),
                                  other.srcrank(), owned, COL_MAJOR),
                      m_parent.context());
        ScalaMat tmp2(LocalMatrix(new double[other.m() * other.n()], other.m(),
                                  other.n(), other.srcrank(), true, COL_MAJOR),
                      m_parent.context());
        return (*this += tmp1.transpose(tmp2).submatrix(
                    rowrange(1, other.m()), colrange(1, other.n())));
    }
}

ScalaMat::Submatrix &ScalaMat::Submatrix::operator=(LocalMatrix &&other)
{
    if (other.ordering() == COL_MAJOR)
    {
        ScalaMat A(std::move(other));
        return (*this = A(rowrange(1, A.m()), colrange(1, A.n())));
    }
    else
    {
        bool owned = other.own_data();
        ScalaMat tmp1(LocalMatrix(other.release_data(), other.n(), other.m(),
                                  other.srcrank(), owned, COL_MAJOR));
        ScalaMat tmp2(LocalMatrix(new double[other.m() * other.n()], other.m(), other.n(),
                                  other.srcrank(), true, COL_MAJOR));
        return (*this = tmp1.transpose(tmp2)(rowrange(1, other.m()), colrange(1, other.n())));
    }
}

ScalaMat::Submatrix::Submatrix(ScalaMat &parent, std::pair<int, int> ibounds,
                               std::pair<int, int> jbounds) : m_parent(parent), m_i(ibounds.first),
                                                              m_j(jbounds.first), m_m(ibounds.second - m_i + 1),
                                                              m_n(jbounds.second - m_j + 1)
{
    if (ibounds.second < ibounds.first || jbounds.second < jbounds.first)
    {
        if (mpi_global_params().first == 0)
        {
            std::cerr << "Error constructing submatrix - second element of "
                      << "index range must be >= the first\n";
            std::cerr << "(in " << __func__ << " at " << __FILE__ << ':'
                      << __LINE__ << ")\n";
        }
        throw ScalaMatException("Submatrix constructor: nonincreasing range");
    }
}

ScalaMat::Submatrix ScalaMat::submatrix(RowRange ibounds, ColRange jbounds)
{
    return Submatrix(*this, std::make_pair(ibounds.lower, ibounds.upper),
                     std::make_pair(jbounds.lower, jbounds.upper));
}

const ScalaMat::Submatrix ScalaMat::submatrix(RowRange ibounds,
                                              ColRange jbounds) const
{
    return Submatrix(const_cast<ScalaMat &>(*this),
                     std::make_pair(ibounds.lower, ibounds.upper),
                     std::make_pair(jbounds.lower, jbounds.upper));
}

LocalMatrix &LocalMatrix::operator=(const ScalaMat::Submatrix &other)
{
    std::array<int, 9> dst_desc, src_desc;
    bool owned = m_own_data;
    m_own_data = false;

    switch (m_ordering)
    {
    case ROW_MAJOR:
    {
        m_ordering = COL_MAJOR;
        std::swap(m_m, m_n);
        bool need_copy = false;
        if (mpi_rank() == m_rank)
        {
            if (!other.parent().context())
            {
                need_copy = true;
            }
        }
        MPI_Bcast(&need_copy, 1, MPI_C_BOOL, m_rank, MPI_COMM_WORLD);

        if (!need_copy)
        {
            ScalaMat target(std::move(*this), other.parent().context());

            target.initialize_descriptor(dst_desc);
            other.parent().initialize_descriptor(src_desc);
            pdgeadd_wrapper('T', target.m(), target.n(), 1.0,
                            other.parent().data(), other.i(), other.j(),
                            src_desc.data(), 0.0, target.data(), 1, 1,
                            dst_desc.data());
        }
        else
        {
            ScalaMat tmp(m_n, m_m);
            tmp(rowrange(1, m_n), colrange(1, m_m)) = other;
            ScalaMat target(std::move(*this));
            tmp.transpose(target);
        }

        std::swap(m_m, m_n);
        m_ordering = ROW_MAJOR;
        break;
    }
    default:
    {
        ScalaMat target(std::move(*this));
        target(rowrange(1, target.m()), colrange(1, target.n())) = other;
        break;
    }
    }

    m_own_data = owned;
    return *this;
}

void ScalaMat::initialize_descriptor(std::array<int, 9> &desc) const noexcept
{
    assert(m_data != nullptr);
    desc[0] = 1;
    desc[1] = m_ctxt ? m_ctxt->id() : -1;
    desc[2] = m_globm;
    desc[3] = m_globn;
    desc[4] = m_mb;
    desc[5] = m_nb;
    desc[6] = m_rowsrc;
    desc[7] = m_colsrc;

    if (m_ctxt)
    {
        int nprow, npcol, pi, pj;
        blacs_gridinfo_wrapper(m_ctxt->id(), &nprow, &npcol, &pi, &pj);
        int mm = numroc_wrapper(m_globm, m_mb, pi, m_rowsrc, nprow);
        desc[8] = std::max(1, mm);
    }
    else
    {
        desc[8] = 1;
    }
}

ScalaMat ScalaMat::Submatrix::similar() const
{
    return ScalaMat(m(), n(), m_parent.m_mb, m_parent.m_nb,
                    m_parent.context(), m_parent.m_rowsrc, m_parent.m_colsrc);
}

void ScalaMat::print_descriptor(std::ostream &os,
                                const std::array<int, 9> &desc)
{
    os << "desc[CTXT_] = " << desc[1] << '\n';
    os << "desc[M_]" << desc[2] << '\n';
    os << "desc[N_]" << desc[3] << '\n';
    os << "desc[MB_]" << desc[4] << '\n';
    os << "desc[NB_]" << desc[5] << '\n';
    os << "desc[ROWSRC_]" << desc[6] << '\n';
    os << "desc[COLSRC_]" << desc[7] << '\n';
    os << "desc[MXLLD_]" << desc[8] << '\n';
}

void ScalaMat::dump_debug_info(std::ostream &os, bool show_data) const
{
    if (!m_ctxt)
    {
        os << "Rank " << mpi_rank() << ": not a member of this matrix's context\n";
        return;
    }

    int nprow, npcol, pi, pj;
    blacs_gridinfo_wrapper(m_ctxt->id(), &nprow, &npcol, &pi, &pj);
    os << "Rank " << mpi_rank() << ": m_ctxt = " << m_ctxt->id()
       << "; " << nprow << " x " << npcol << " process grid; I am process ("
       << pi << ", " << pj << ")\n";

    os << "Blocking factors: (" << m_mb << ", " << m_nb << "); "
       << "matrix is distributed starting from (" << m_rowsrc << ", "
       << m_colsrc << ")\n";
    int mm = numroc_wrapper(m_globm, m_mb, pi, m_rowsrc, nprow);
    int mn = numroc_wrapper(m_globn, m_nb, pj, m_colsrc, npcol);
    os << "My local data array is " << mm << " x " << mn;

    if (show_data && mm > 0 && mn > 0)
    {
        auto prec = os.precision();
        os.precision(4);
        os << ":\n";

        for (int i = 0; i < mm; ++i)
        {
            for (int j = 0; j < mn; ++j)
            {
                os << m_data[j * mm + i] << ' ';
            }
            os << '\n';
        }

        os.precision(prec);
    }
    else
    {
        os << ".\n";
    }
}

} /* namespace ScalaWRAP */
