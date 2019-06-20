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
       * @pre dim > 0
       * @pre samples_per_time_interval > 0
       *
       * @param[in] dim The dimension of the system on this processor.
       * @param[in] samples_per_time_interval The maximum number of samples in
       *                                      each time interval.
       * @param[in] basis_file_name The base part of the name of the file
       *                            containing the basis vectors.  Each process
       *                            will append its process ID to this base
       *                            name.
       * @param[in] max_basis_dimension (typemax(int)) The maximum number of
       *                                vectors returned in the basis.
       * @param[in] sigma_tolerance This tolerance is based on the ratio of
       *                            singular values to the largest singular
       *                            value. If sigma[i] / sigma[0] < sigma_tolerance,
       *                            the associated vector is dropped from the
       *                            basis.
       * @param[in] debug_algorithm If true results of static svd algorithm
       *                            will be printed to facilitate debugging.
       * @param[in] file_format The format of the file containing the basis
       *                        vectors.
       */
      StaticSVDBasisGenerator(
         int dim,
         int samples_per_time_interval,
         const std::string& basis_file_name,
         int max_basis_dimension = std::numeric_limits<int>::max(),
         double sigma_tolerance = 0,
         bool debug_algorithm = false,
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
