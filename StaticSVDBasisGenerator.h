/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The concrete wrapper class for static SVD algorithm and
//              sampler.  Implements interface of SVDBasisGenerator.

#ifndef included_StaticSVDBasisGenerator_h
#define included_StaticSVDBasisGenerator_h

#include "SVDBasisGenerator.h"
#include "StaticSVD.h"

#include <limits>

namespace CAROM {

/**
 * Class StaticSVDBasisGenerator implements the interface of base class
 * SVDBasisGenerator for the static svd algorithm.
 */
class StaticSVDBasisGenerator : public SVDBasisGenerator
{
   public:
      /**
       * @brief Constructor.
       *
       * @param[in] options The struct containing the options for this basis
       *                    generator.
       * @param[in] basis_file_name The base part of the name of the file
       *                            containing the basis vectors.  Each process
       *                            will append its process ID to this base
       *                            name.
       * @param[in] file_format The format of the file containing the basis
       *                        vectors.
       */
      StaticSVDBasisGenerator(
         StaticSVDOptions options,
         const std::string& basis_file_name = "",
         Database::formats file_format = Database::HDF5);

      /**
       * @brief Destructor.
       */
      virtual
      ~StaticSVDBasisGenerator();

      bool
      updateRightSV() { return true; }

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      StaticSVDBasisGenerator();

      /**
       * @brief Unimplemented copy constructor.
       */
      StaticSVDBasisGenerator(
         const StaticSVDBasisGenerator& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      StaticSVDBasisGenerator&
      operator = (
         const StaticSVDBasisGenerator& rhs);
};

}

#endif
