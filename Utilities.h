/******************************************************************************
 *
 * Copyright (c) 2013-2016, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 * Written by William Arrighi wjarrighi@llnl.gov
 * CODE-686965
 * All rights reserved.
 *
 * This file is part of libROM.
 * For details, see https://computation.llnl.gov/librom
 * Please also read README_BSD_NOTICE.
 *
 * Redistribution and use in source and binary forms, with or without
 * modifications, are permitted provided that the following conditions are met:
 *
 *    o Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the disclaimer below.
 *    o Redistribution in binary form must reproduce the above copyright
 *      notice, this list of conditions and the disclaimer (as noted below) in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *    o Neither the name of the LLNS/LLNL nor the names of its contributors may
 *      be used to endorse or promote products derived from this software
 *      without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
 * LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OR SUCH DAMAGE.
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
 * unused variables.
 */
#define CAROM_NULL_USE(variable)                         \
   do {                                                  \
      if (0) { char* temp = (char *)&variable; temp++; } \
   } while (0)

/**
 * @brief Throw an error assertion from within any C++ source code.
 *
 * The macro argument may be any standard ostream expression.  The file and
 * line number of the abort are also printed.
 */
#define CAROM_ERROR(X)                                       \
   do {                                                      \
      std::ostringstream os;                                 \
      os << X << std::ends;                                  \
      CAROM::Utilities::abort(os.str(), __FILE__, __LINE__); \
   } while (0)

/**
 * @brief Throw an error assertion from within any C++ source code if the given
 * expression is not true.
 *
 * This is a parallel-friendly version of assert.  The file and line number of
 * the abort are also printed.
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
 * Utilities contains basic, static, utility routines for error reporting,
 * string manipulations, etc.
 */
struct Utilities
{
   /**
    * @brief Cleanly ends the program when something horrible happend and
    * prints a message about what took place.
    *
    * Takes into account whether MPI is or isn't running to decide how to die.
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
    * @brief Converts a processor ID to a string.
    *
    * Use this to ensure same width is used when converting a processor ID to
    * a string representation.
    *
    * @param[in] processorID of the processor
    *
    * @return The string representation of processor ID of fixed width
    * prepended with 0s.
    */
   static std::string
   processorToString(
      int processorID);
};

}

#endif
