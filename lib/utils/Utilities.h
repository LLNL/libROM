/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Utility functions for error reporting, etc.

#ifndef included_utilities_h
#define included_utilities_h

#include "CAROM_config.h"
#include <iostream>
#include <sstream>
#include <string>

namespace CAROM {

/**
 * @brief A null use of a variable, use to avoid GNU compiler warnings about
 *        unused variables.
 */
#define CAROM_NULL_USE(variable)                         \
   do {                                                  \
      if (0) { char* temp = (char *)&variable; temp++; } \
   } while (0)

/**
 * @brief Throw an error assertion from within any C++ source code.
 *        The macro argument may be any standard ostream expression. The file
 *        and line number of the abort are also printed.
 */
#define CAROM_ERROR(X)                                       \
   do {                                                      \
      std::ostringstream os;                                 \
      os << X << std::ends;                                  \
      CAROM::Utilities::abort(os.str(), __FILE__, __LINE__); \
   } while (0)

/**
 * @brief Throw an error assertion from within any C++ source code if the given
 *        expression is not true. This is a parallel-friendly version of assert.
 *        The file and line number of the abort are also printed.
 */
#define CAROM_VERIFY(EXP)                                       \
    do {                                                         \
       if (!(EXP)) {                                             \
          std::ostringstream os;                                 \
          os << "Failed verify: " << # EXP << std::ends;      \
          CAROM::Utilities::abort(os.str(), __FILE__, __LINE__); \
       }                                                         \
    } while (0)

/**
 * @brief Throw an error assertion from within any C++ source code if the given
 *        expression is not true. This is a parallel-friendly version of assert.
 *        The file and line number of the abort are also printed.
 */
#ifdef DEBUG_CHECK_ASSERTIONS

#define CAROM_ASSERT(EXP)                                       \
   do {                                                         \
      if (!(EXP)) {                                             \
         std::ostringstream os;                                 \
         os << "Failed assertion: " << # EXP << std::ends;      \
         CAROM::Utilities::abort(os.str(), __FILE__, __LINE__); \
      }                                                         \
   } while (0)
#else

/*
 * No assertion checking
 */
#define CAROM_ASSERT(EXP)

#endif

/**
 * Struct BasisGenerator defines Utilities contains basic, static, utility
 * routines for error reporting, string manipulations, etc.
 */
struct Utilities
{
    /**
     * @brief Cleanly ends the program when something horrible happened and
     *        prints a message about what took place. Takes into account
     *        whether MPI is or isn't running to decide how to die.
     *
     * @param[in] message Message to print about the cause of the abort.
     * @param[in] filename Name of the file where the abort was called.
     * @param[in] line Line number in the file where the abort was called.
     */
    static void
    abort(
        const std::string& message,
        const std::string& filename,
        int line);

    /**
     * @brief Converts a processor ID to a string. Use this to ensure same width
     *        is used when converting a processor ID to a string representation.
     *
     * @param[in] processorID of the processor
     *
     * @return The string representation of processor ID of fixed width
     * prepended with 0s.
     */
    static std::string
    processorToString(
        int processorID);

    /**
     * @brief Verifies if a file exists.
     *
     * @param[in] filename Name of the file to be verified.
     *
     * @return true if the file exists.
     */
    static bool
    file_exist(
        const std::string& filename);
};

}

#endif
