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
<<<<<<< HEAD
       * @param[in] save_state If true the state of the SVD will be written to
       *                       disk when the object is deleted.  If there are
       *                       multiple time intervals then the state will not
       *                       be saved as restoring such a state makes no
       *                       sense.
       * @param[in] restore_state If true the state of the SVD will be restored
       *                          when the object is created.
       * @param[in] update_rightSV Whether to update the right SV or not.
       * @param[in] write_snapshots Whether to automatically write snapshots matrices
       *                        instead of basis matrices.
=======
>>>>>>> origin/master
       * @param[in] file_format The format of the file containing the basis
       *                        vectors.
       * @param[in] write_snapshots Whether to automatically write snapshots matrices
       *                        instead of basis matrices.
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
