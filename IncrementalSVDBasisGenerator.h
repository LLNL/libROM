/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The concrete wrapper class for a specific incremental SVD
//              algorithm and sampler.  Implements interface of
//              SVDBasisGenerator.

#ifndef included_IncrementalSVDBasisGenerator_h
#define included_IncrementalSVDBasisGenerator_h

#include "SVDBasisGenerator.h"
#include "IncrementalSVD.h"

namespace CAROM {

/**
 * Class IncrementalSVDBasisGenerator implements the interface of base class
 * SVDBasisGenerator for the incremental svd algorithm.  Either the fast update
 * or the standard incremental algorithm may be specified through the
 * constructor.
 */

class IncrementalSVDBasisGenerator : public SVDBasisGenerator
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
      IncrementalSVDBasisGenerator(
         IncrementalSVDOptions options,
         const std::string& basis_file_name = "",
         Database::formats file_format = Database::HDF5);

      /**
       * @brief Destructor.
       */
      virtual
      ~IncrementalSVDBasisGenerator();

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      IncrementalSVDBasisGenerator();

      /**
       * @brief Unimplemented copy constructor.
       */
      IncrementalSVDBasisGenerator(
         const IncrementalSVDBasisGenerator& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      IncrementalSVDBasisGenerator&
      operator = (
         const IncrementalSVDBasisGenerator& rhs);
};

}

#endif
