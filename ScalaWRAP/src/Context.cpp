/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

#include "Context.hpp"
#include "C_interface.h"
#include "MPI_utils.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <set>
#include <tuple>
#include <utility>

namespace ScalaWRAP {

struct ContextException {
    const char* msg;
    const char* what() const noexcept
    {
        return msg;
    }

    ContextException(const char* m) : msg(m) {}
};

namespace {

std::shared_ptr<const Context> DEFAULT_CONTEXT_PTR = nullptr;

std::pair<int, int> best_distribution(int nprocs)
{
    int i = std::floor(std::sqrt(nprocs));
    while (nprocs % i != 0) {
        i -= 1;
    }
    return std::make_pair(i, nprocs / i);
}

}

Context::Context(int nprow, int npcol)
{
    int nprocs, mrank;
    std::tie(nprocs, mrank) = mpi_global_params();
    if (nprocs < nprow * npcol)
    {
        if (mrank == 0)
        {
            std::cerr << "Error: request for " << nprow << " x " << npcol
                      << " process grid, only " << nprocs << " processes available\n";
            std::cerr << "(in " << __func__ << " at " << __FILE__ << ':' 
                      << __LINE__ << '\n';
        }
        throw ContextException("Context constructor: Requested unfulfillable grid");
    }
    else
    {
        m_contextid = sl_init_wrapper(nprow, npcol);
    }
}

Context::Context(const std::vector<int>& process_map, int nprow)
{
    int mrank, nprocs;
    std::tie(nprocs, mrank) = mpi_global_params();
    
    assert(nprocs >= 1);
    if (process_map.size() > static_cast<unsigned>(nprocs)) {
        if (mrank == 0) {
            std::cerr << "Error: process_map.size() = " << process_map.size()
                    << " but size of MPI_COMM_WORLD is only " << nprocs << '\n';
            std::cerr << "(in " << __func__ << " at " << __FILE__ << ':'
                    << __LINE__ << '\n';
        }
        throw ContextException("Context process_map constructor: Requested unfulfillable grid");
    }

    // Check for multiple ranks in process_map, or processes outside of the
    // valid index range.
    std::set<int> rankset;
    for (auto rank: process_map) {
        if (rankset.count(rank) == 1) {
            if (mrank == 0) {
                std::cerr <<
                    "Error: Multiple entries in process_map refer to the same rank\n"
                    << "(in " << __func__ << " at " << __FILE__ << ':'
                    << __LINE__ << '\n';
            }
            throw ContextException("Context process_map constructor: multiple specs of same rank");
        } else if (rank >= nprocs) {
            if (mrank == 0) {
                std::cerr <<
                    "Error: rank = " << rank << " is outside of the index "
                    << "range for processes (nprocs = " << nprocs << ")\n";
                std::cerr << "(in " << __func__ << " at " << __FILE__ << ':'
                    << __LINE__ << '\n';
            }
            throw ContextException("Context process_map constructor: rank outside the valid range");
        }
        rankset.insert(rank);
    }

    if (process_map.size() % nprow != 0) {
        if (mrank == 0) {
            std::cerr << 
                "Error: process_map.size() is not a multiple of nprow\n"
                << "(in " << __func__ << " at " << __FILE__ << ':' << __LINE__
                << '\n';
        }
        throw ContextException("Context process_map constructor: process_map.size % nprow != 0");
    }

    m_contextid = blacs_gridmap_wrapper(process_map.data(), nprow, nprow,
                                        process_map.size() / nprow);
}

std::shared_ptr<const Context> Context::default_context()
{
    if (!DEFAULT_CONTEXT_PTR) {
        int nprow, npcol;
        std::tie(nprow, npcol) = best_distribution(mpi_global_params().first);
        DEFAULT_CONTEXT_PTR = std::make_shared<Context>(nprow, npcol);
    }
    return DEFAULT_CONTEXT_PTR;
}

void finalize_scalawrap()
{
    DEFAULT_CONTEXT_PTR && (DEFAULT_CONTEXT_PTR = nullptr);
}

std::tuple<int, int, int, int> Context::getinfo() const
{
    int nprow, npcol, pi, pj;
    blacs_gridinfo_wrapper(m_contextid, &nprow, &npcol, &pi, &pj);
    return std::make_tuple(nprow, npcol, pi, pj);
}

} /* namespace ScalaWRAP */