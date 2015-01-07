/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: The concrete wrapper class for a specific incremental SVD
 *              algorithm and sampler.  Implements interface of
 *              SVDBasisGenerator.
 *
 *****************************************************************************/

#ifndef included_IncrementalSVDBasisGenerator_h
#define included_IncrementalSVDBasisGenerator_h

#include "SVDBasisGenerator.h"

namespace CAROM {

/**
 * Class IncrementalSVDBasisGenerator implements the interface of base class
 * SVDBasisGenerator for the incremental svd algorithm.  Either the fast update
 * or the naive incremental algorithm may be specified through the constructor.
 */
class IncrementalSVDBasisGenerator : public SVDBasisGenerator
{
   public:
      /**
       * @brief Constructor.
       *
       * @pre dim > 0
       * @pre redundancy_tol > 0.0
       * @pre samples_per_time_interval > 0
       * @pre sampling_tol > 0.0
       * @pre max_time_between_samples > 0.0
       *
       * @param[in] dim The dimension of the system on this processor.
       * @param[in] redundancy_tol Tolerance to determine if a sample is
       *                           redundant or not.
       * @param[in] skip_redundant If true skip redundant samples.
       * @param[in] samples_per_time_interval The maximum number of samples in
       *                                      each time interval.
       * @param[in] sampling_tol Sampling control tolerance.  Limits error in
       *                         projection of solution into reduced order
       *                         space.
       * @param[in] max_time_between_samples Hard upper bound on sampl.
       * @param[in] fast_update If true use the fast update algorithm.
       * @param[in] basis_file_name The base part of the name of the file
       *                            containing the basis vectors.  Each process
       *                            will append its process ID to this base
       *                            name.
       * @param[in] debug_rom If true results of incremental svd algorithm
       *                      will be printed to facilitate debugging.
       * @param[in] file_format The format of the file containing the basis
       *                        vectors.
       */
      IncrementalSVDBasisGenerator(
         int dim,
         double redundancy_tol,
         bool skip_redundant,
         int samples_per_time_interval,
         double sampling_tol,
         double max_time_between_samples,
         bool fast_update,
         const std::string& basis_file_name,
         bool debug_rom = false,
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
