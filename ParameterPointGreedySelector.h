/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: This class greedily selects parameter points
//              for the construction of a ROM database.

#ifndef included_ParameterPointGreedySelector_h
#define included_ParameterPointGreedySelector_h

#include "BasisGenerator.h"
#include <algorithm>

namespace CAROM {

class BasisGenerator;
class Matrix;
class Vector;

/**
 * Class ParameterPointGreedySelector is an class defining the interface for
 * the generation of basis vectors via the svd method.  This class wraps the
 * SVD algorithm and sampler and controls all aspects
 * of basis vector generation.
 */
class ParameterPointGreedySelector
{
   public:
     /**
      * @brief Constructor.
      *
      * @param[in] options The struct containing the options for this basis
      *                    generator.
      * @param[in] incremental Whether to conduct static or incremental SVD
      * @param[in] basis_file_name The base part of the name of the file
      *                            containing the basis vectors.  Each process
      *                            will append its process ID to this base
      *                            name.
      * @param[in] file_format The format of the file containing the basis
      *                        vectors.
      */
     ParameterPointGreedySelector(
        std::vector<Vector> parameter_points,
        double tolerance,
        double saturation,
        int subset_size,
        int convergence_subset_size,
        int random_seed);

    ParameterPointGreedySelector(
        std::vector<double> parameter_points,
        double tolerance,
        double saturation,
        int subset_size,
        int convergence_subset_size,
        int random_seed);

      /**
       * @brief Destructor.
       */
      ~ParameterPointGreedySelector();

      /**
       * @brief Returns the next parameter point for which sampling is required.
       *
       * @return The index of the next parameter point, or -1 if the greedy
       *         procedure has terminated.
       */
      int
      computeNextSampleParameterPoint();

      /**
       * @brief Returns the nearest local ROM to the specified parameter point.
       *
       * @return A pointer to the BasisGenerator containing the local ROM
       */
      const BasisGenerator*
      getNearestLocalROM(double parameter_point);

      /**
       * @brief Add rom to database.
       *
       * @param[in] rom A pointer to the rom (not owned by ParameterPointGreedySelector).
       *
       */
      void
      addROMToDatabase(BasisGenerator* rom);

  private:

      /**
       * @brief The parameter points to explore.
       */
      std::vector<Vector> d_parameter_points;

      /**
       * @brief The parameter point indices (used to generate the random subsets)
       */
      std::vector<int> d_parameter_point_random_indices;

      /**
       * @brief The current errors of the parameter points
       */
      std::vector<double> d_parameter_point_errors;

      /**
       * @brief The current max error of the parameter points
       */
      double d_max_error;

     /**
      * @brief The convergence tolerance used to terminate the greedy procedure.
      */
      double d_tol;

     /**
      * @brief The saturation constant
      */
      double d_sat;

     /**
      * @brief The size of the subset of parameter points used per iteration.
      */
      int d_subset_size;

     /**
      * @brief The size of the subset of parameter points used to check convergence.
      */
      int d_convergence_subset_size;


    /**
     * @brief The next point to add to the database.
     */
     int d_next_point;

      /**
       * @brief Whether the database is ready to accept a ROM for the current
       *        iteration.
       */
       bool d_iteration_completed;

       /**
        * @brief Whether the greedy procedure has completed.
        */
        bool d_procedure_completed;

        /**
         * @brief Random engine used to generate subsets
         */
       std::default_random_engine rng;

      /**
       * @brief The ROM database.
       */
      std::vector<BasisGenerator*> d_rom_database;
};

}

#endif
