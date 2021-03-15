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
//              for the construction of a ROM database. The implemented
//              greedy algorithm is obtained from algorithm 2 of Choi et. al's
//              paper "Gradient-based constrained optimization using a database
//              of linear reduced-order models".
//
//              The greedy algorithm workflow is as follows:
//              1. Construct the GreedyParameterPointSelector by giving it a
//                 domain of parameter points.
//              2. Request a parameter point to sample.
//              3. Request a parameter point to compute an error residual for.
//              4. Request the nearest ROM to the parameter point requiring a
//                 residual.
//              5. Give the computed residual to the GreedyParameterPointSelector.
//              6. Repeat steps 4 and 5 until the GreedyParameterPointSelector
//                 no longer requires more residuals to be computed.
//              7. Repeat steps 2 to 6 until the GreedyParameterPointSelector
//                 no longer requires more parameter points to be sampled.
//              8. The ROM database is now complete, meeting the error tolerance
//                 for all parameter points within the domain.

#ifndef included_GreedyParameterPointSelector_h
#define included_GreedyParameterPointSelector_h

#include "BasisGenerator.h"
#include <random>
#include <set>
#include <algorithm>

namespace CAROM {

class BasisGenerator;
class Matrix;
class Vector;

/**
 * Class GreedyParameterPointSelector is an class defining the interface a
 *       greedy algorithm that given a domain of parameter points, iteratively
 *       returns the next best parameter point to sample in order to create
 *       a ROM database efficiently.
 */
class GreedyParameterPointSelector
{
   public:
     /**
      * @brief Constructor.
      *
      * @param[in] parameter_points A vector of CAROM::Vectors containing
                                    the different parameter points.
      * @param[in] tolerance A tolerance value for which to end the algorithm.
      * @param[in] saturation A saturation constant.
      * @param[in] subset_size The size of the random subset.
      * @param[in] convergence_subset_size The size of the convergence subset.
      * @param[in] use_centroid Whether to use the centroid heuristic when
                                determining the first parameter point to sample.
      * @param[in] random_seed A random seed.
      * @param[in] debug_algorithm Whether to turn off all randomness for
      *                            debugging purposes.
      */
     GreedyParameterPointSelector(
        std::vector<Vector> parameter_points,
        double tolerance,
        double saturation,
        int subset_size,
        int convergence_subset_size,
        bool use_centroid = true,
        int random_seed = 1,
        bool debug_algorithm = false);

    GreedyParameterPointSelector(
        std::vector<double> parameter_points,
        double tolerance,
        double saturation,
        int subset_size,
        int convergence_subset_size,
        bool use_centroid = true,
        int random_seed = 1,
        bool debug_algorithm = false);

    GreedyParameterPointSelector(
        double param_space_min,
        double param_space_max,
        double param_space_size,
        double tolerance,
        double saturation,
        int subset_size,
        int convergence_subset_size,
        bool use_centroid = true,
        int random_seed = 1,
        bool debug_algorithm = false);

      /**
       * @brief Destructor.
       */
      ~GreedyParameterPointSelector();

      /**
       * @brief Returns the next parameter point for which sampling is required.
       *
       * @return The index of the point in the list of parameters, or -1 if a point
       *         is not ready to be sampled.
       */
      int
      getNextParameterPoint();

      /**
       * @brief Compute point whose residual the greedy algorithm requires.
       *
       * @param[in]
       *
       * @return The index of the point in the list of parameters, or -1 if a
       *         residual is not currently required.
       */
      int
      getNextPointRequiringResidual();

      /**
       * @brief Set the residual error of the specified parameter point.
       *
       * @return The index of the point in the list of parameters.
       */
      void
      setPointResidual(double error, int rank, int num_procs);

      /**
       * @brief Returns the index to the nearest local ROM to the specified parameter point.
       *
       * @return The index of the point in the list of parameters.
       */
      int
      getNearestROM(int index);

      /**
       * @brief Get the domain of the parameter points.
       */
      std::vector<Vector>
      getParameterPointDomain();

      /**
       * @brief Print the sampled parameter points to a file.
       */
      void
      printSampledPoints(std::string const& path);

  private:

      void constructObject(
          std::vector<Vector> parameter_points,
          double tolerance,
          double saturation,
          int subset_size,
          int convergence_subset_size,
          bool use_centroid,
          int random_seed,
          bool debug_algorithm);

      /**
       * @brief The parameter points to explore.
       */
      std::vector<Vector> d_parameter_points;

      /**
       * @brief The parameter points that were already sampled.
       */
      std::set<int> d_parameter_sampled_indices;

      /**
       * @brief The parameter point indices (used to generate the random subsets)
       */
      std::vector<int> d_parameter_point_random_indices;

      /**
       * @brief The current errors of the parameter points
       */
      std::vector<double> d_parameter_point_errors;

      /**
       * @brief The current max error of the parameter points of the current iteration
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
      * @brief The next parameter point to sample.
      */
      int d_next_point_to_sample;

     /**
      * @brief The next parameter point requiring a residual.
      */
      int d_next_point_requiring_residual;

     /**
      * @brief Whether the use the centroid heuristic for obtaining the first
      *        parameter point.
      */
      bool d_use_centroid;

      /**
       * @brief Whether the database has already computed a new parameter point
       *        for the current iteration.
       */
       bool d_iteration_started;

       /**
        * @brief Whether the database has already computed a new paramter point
        *        requiring a residual.
        */
        bool d_point_requiring_residual_computed;

       /**
        * @brief Whether the database has already created a random subset
        *        for this iteration.
        */
        bool d_subset_created;

        /**
         * @brief Turn off randomness for debugging purposes.
         *
         */
        bool d_debug_algorithm;

        /**
         * @brief An internal counter.
         *
         */
       int d_counter;

       /**
        * @brief An internal subset counter.
        *
        */
       int d_subset_counter;

       /**
        * @brief Whether the greedy procedure has completed.
        */
        bool d_procedure_completed;

        /**
         * @brief Random engine used to generate subsets
         */
       std::default_random_engine rng;
};

// Given a a vector/double, find the nearest point in a domain.
int getNearestPoint(std::vector<Vector> paramPoints, Vector point);
int getNearestPoint(std::vector<double> paramPoints, double point);

}

#endif
