/*******************************************************************************
 *
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC under the terms
 * of the MIT license. See the top-level COPYRIGHT file for details.
 * 
 * SPDX-License-Identifier: MIT
 * 
 ******************************************************************************/

/*!
 * @file Context.hpp
 * 
 * @brief This is the public interface to the `Context` wrapper class used to
 * manage BLACS contexts in this interface. The `ScalaMAT` class depends
 * directly on this interface.
 */

#ifndef CONTEXT_HPP
#define CONTEXT_HPP

#include "C_interface.h"

#include <memory>
#include <tuple>
#include <vector>

namespace ScalaWRAP {

/*!
 * @brief Class providing encapsulation for the raw integer handle that BLACS
 * uses to manage "contexts"
 * 
 * The BLACS concept of a "context" in which a distributed matrix lives is
 * fundamental to writing code using ScaLAPACK. The context is represented by
 * an opaque integer handle. This class exists to protect code handling
 * ScaLAPACK matrices from accidentally deleting or changing a context and to
 * provide a nice interface to the BLACS functions that manipulate contexts.
 * 
 * It is not possible for the user to directly construct a `Context`. Rather,
 * you should always handle a `shared_ptr<Context>`. This is because when the
 * `Context` destructor is called, if there is an associated BLACS context it
 * will be freed. Using a `shared_ptr` to a single instance of `Context` for
 * all matrices sharing a context makes sure that the context is not released
 * prematurely. Please don't try to circumvent this interface by dereferencing
 * the smart pointers returned by the interface, although of course it is
 * possible.
 */
class Context
{
public:
    /*! @brief Retrieve the integer handle used by BLACS */
    int id() const noexcept { return m_contextid; }

    /*! @brief If this object has an associated BLACS context, free it. */
    ~Context()
    {
        if (m_contextid != -1) {
            blacs_gridexit_wrapper(m_contextid);
        }
    }

    /*!
     * @brief Returns a default global context.
     * 
     * The default context lives from the first call to this function until the
     * program exits. It maps all of the processes in MPI_COMM_WORLD into a
     * process grid whose size is as close to square as possible. When the grid
     * must be rectangular, more processes are put in the row direction.
     */
    static std::shared_ptr<const Context> default_context();

    /*!
     * @brief Construct a context with an `nprow` by `npcol` grid.
     */
    static std::shared_ptr<const Context> make_context(int nprow, int npcol)
    {
        return std::make_shared<Context>(nprow, npcol);
    }

    /*!
     * @brief Construct a context with a manually specified mapping of
     * processes to the grid.
     * 
     * @param[in] process_map A vector containing the ranks of processes to be
     * mapped to the grid. A process may not appear twice in this specification,
     * and the length of the vector must be a multiple of `nprow`. The data is
     * interpreted as a column major array, as the raw pointer is passed
     * directly to a Fortran wrapper routine.
     * 
     * @param[in] nprow The number of rows in the created process grid. The
     * number of columns is `process_map.size() / nprow`
     * 
     * @note Throws an exception if an illegal `process_map` is given.
     */
    static std::shared_ptr<const Context> make_context(
        const std::vector<int>& process_map, int nprow)
    {
        return std::make_shared<Context>(process_map, nprow);
    }

    Context(int nprow, int npcol);
    Context(const std::vector<int>& process_map, int nprow);

    /*!
     * @brief Get a tuple `(nprow, npcol, pi, pj)` of the parameters of the
     * active context.
     * 
     * This operation corresponds to a call to blacs_gridinfo. `pi` and `pj` are
     * the zero-based indices of this process's position in the process grid.
     */
    std::tuple<int, int, int, int> getinfo() const;

private:
    int m_contextid;
};

void finalize_scalawrap();

} /* namespace ScalaWRAP */

#endif /* CONTEXT_HPP */
