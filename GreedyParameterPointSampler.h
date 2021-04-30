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
//              of linear reduced-order models": https://arxiv.org/abs/1506.07849
//
//              The greedy algorithm workflow is as follows:
//              1. Construct the GreedyParameterPointSampler by giving it a
//                 domain of parameter points.
//              2. Request a parameter point to sample.
//              3. Request a parameter point to compute an error residual for.
//              4. Request the nearest ROM to the parameter point requiring a
//                 residual.
//              5. Give the computed residual to the GreedyParameterPointSampler.
//              6. Repeat steps 4 and 5 until the GreedyParameterPointSampler
//                 no longer requires more residuals to be computed.
//              7. Repeat steps 2 to 6 until the GreedyParameterPointSampler
//                 no longer requires more parameter points to be sampled.
//              8. The ROM database is now complete, meeting the error tolerance
//                 for all parameter points within the domain.

#ifndef included_GreedyParameterPointSampler_h
#define included_GreedyParameterPointSampler_h

#include "Vector.h"
#include <string>
#include <vector>
#include <set>
#include <random>

/* Use C++11 built-in shared pointers if available; else fallback to Boost. */
#if __cplusplus >= 201103L
#include <memory>
#else
#include <boost/shared_ptr.hpp>
#endif

namespace CAROM {

class Vector;

struct GreedyResidualPoint {
    std::shared_ptr<Vector> point;
    std::shared_ptr<Vector> localROM;
};

/**
 * Class GreedyParameterPointSampler is a class defining the interface of a
 *       greedy algorithm that given a domain of parameter points, iteratively
 *       returns the next best parameter point to sample in order to create
 *       a ROM database efficiently.
 */
class GreedyParameterPointSampler
{
public:
    /**
     * @brief Constructor.
     *
     * @param[in] parameter_points A vector of CAROM::Vectors containing
                                   the different parameter points.
     * @param[in] check_local_rom Compute local ROM residual each iteration.
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
    GreedyParameterPointSampler(
        std::vector<Vector> parameter_points,
        bool check_local_rom,
        double relative_error_tolerance,
        double saturation,
        int subset_size,
        int convergence_subset_size,
        std::string output_log_path = "",
        std::string warm_start_file_name = "",
        bool use_centroid = true,
        int random_seed = 1,
        bool debug_algorithm = false);

    GreedyParameterPointSampler(
        std::vector<double> parameter_points,
        bool check_local_rom,
        double relative_error_tolerance,
        double saturation,
        int subset_size,
        int convergence_subset_size,
        std::string output_log_path = "",
        std::string warm_start_file_name = "",
        bool use_centroid = true,
        int random_seed = 1,
        bool debug_algorithm = false);

    GreedyParameterPointSampler(
        Vector param_space_min,
        Vector param_space_max,
        int num_parameter_points,
        bool check_local_rom,
        double relative_error_tolerance,
        double saturation,
        int subset_size,
        int convergence_subset_size,
        std::string output_log_path = "",
        std::string warm_start_file_name = "",
        bool use_centroid = true,
        int random_seed = 1,
        bool debug_algorithm = false);

    GreedyParameterPointSampler(
        double param_space_min,
        double param_space_max,
        int num_parameter_points,
        bool check_local_rom,
        double relative_error_tolerance,
        double saturation,
        int subset_size,
        int convergence_subset_size,
        std::string output_log_path = "",
        std::string warm_start_file_name = "",
        bool use_centroid = true,
        int random_seed = 1,
        bool debug_algorithm = false);

    GreedyParameterPointSampler(
        std::string base_file_name,
        std::string output_log_path = "");

    /**
     * @brief Destructor.
     */
    ~GreedyParameterPointSampler();

    /**
     * @brief Returns the next parameter point for which sampling is required.
     *
     * @return The point, or NULL if a point is not ready to be sampled.
     */
    std::shared_ptr<Vector>
    getNextParameterPoint();

    struct GreedyResidualPoint
    getNextPointRequiringRelativeError();

    /**
     * @brief Returns the next parameter point whose residual the greedy
     *        algorithm requires.
     *
     * @param[in]
     *
     * @return A tuple where the first element is the point, and the second
     *         element is the nearest local ROM. Both elements will be NULL if a
     *         residual is not currently required.
     */
    struct GreedyResidualPoint
    getNextPointRequiringResidual();

    void
    setPointRelativeError(double error);

    /**
     * @brief Set the residual error of the specified parameter point.
     */
    void
    setPointResidual(double error, int vec_size);

    int
    getNearestNonSampledPoint(Vector point);

    /**
     * @brief Returns the nearest local ROM to the specified parameter point.
     *
     * @return The point in the list of parameters.
     */
    std::shared_ptr<Vector>
    getNearestROM(Vector point);

    /**
     * @brief Get the domain of the parameter points.
     */
    std::vector<Vector>
    getParameterPointDomain();

    /**
     * @brief Get the sampled parameter points.
     */
    std::vector<Vector>
    getSampledParameterPoints();

    /**
     * @brief Save the object state to a file.
     */
    virtual void
    save(std::string base_file_name);

    /**
     * @brief Check if the greedy algorithm procedure is complete.
     */
    bool
    isComplete();

protected:

    void addDatabaseFromFile(
        std::string const& warm_start_file_name);

    virtual void load(
        std::string base_file_name);

    virtual void constructParameterPoints();

    void constructObject(
        bool check_local_rom,
        double relative_error_tolerance,
        double saturation,
        int subset_size,
        int convergence_subset_size,
        std::string output_log_path,
        bool use_centroid,
        int random_seed,
        bool debug_algorithm);

    void initializeParameterPoints();

    struct GreedyResidualPoint getNextSubsetPointRequiringResidual();

    struct GreedyResidualPoint getNextConvergencePointRequiringResidual();

    std::vector<Vector> generateRandomPoints(int num_points);

    void printResidual(Vector residualPoint, double proc_errors);

    void printErrorIndicatorToleranceNotMet();

    void setSubsetResidual(double proc_errors);

    void setConvergenceResidual(double proc_errors);

    void generateConvergenceSubset();

    void startConvergence();

    virtual void getNextParameterPointAfterConvergenceFailure() = 0;

    /**
     * @brief Returns the index to the nearest local ROM to the specified parameter point.
     *
     * @return The index of the point in the list of parameters.
     */
    int
    getNearestROMIndex(int index, bool ignore_self);

    /**
     * @brief The parameter points to explore.
     */
    std::vector<Vector> d_parameter_points;

    /**
     * @brief The convergence parameter points to explore.
     */
    std::vector<Vector> d_convergence_points;

    /**
     * @brief The minimum of the parameter space.
     */
    Vector d_min_param_point;

    /**
     * @brief The maximum of the parameter space.
     */
    Vector d_max_param_point;

    /**
     * @brief The parameter points that were already sampled.
     */
    std::set<int> d_parameter_sampled_indices;

    /**
     * @brief The parameter point indices (used to generate the random subsets).
     */
    std::vector<int> d_parameter_point_random_indices;

    /**
     * @brief The current errors of the parameter points.
     */
    std::vector<double> d_parameter_point_errors;

    /**
     * @brief The local ROMs used to obtain the current errors of the parameter points.
     */
    std::vector<int> d_parameter_point_local_rom;

    /**
     * @brief Output log path.
     */
    std::string d_output_log_path;

    /**
     * @brief The current max error of the parameter points of the current iteration.
     */
    double d_max_error;

    /**
     * @brief The convergence tolerance used to terminate the greedy procedure.
     */
    double d_curr_relative_error;

    /**
     * @brief The convergence tolerance used to terminate the greedy procedure.
     */
    double d_alpha;

    /**
     * @brief The convergence tolerance used to terminate the greedy procedure.
     */
    double d_error_indicator_tol;

    /**
     * @brief The convergence tolerance used to terminate the greedy procedure.
     */
    double d_relative_error_tol;

    /**
     * @brief The maximum number of parameter points.
     */
    int d_num_parameter_points;

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
     * @brief Whether the check the last sampled local ROM's residual
     *        after each iteration.
     */
    bool d_check_local_rom;

    /**
     * @brief Whether the database has already computed a new parameter point
     *        for the current iteration.
     */
    bool d_iteration_started;

    /**
     * @brief Whether the database is in the convergence verifying phase.
     */
    bool d_convergence_started;

    /**
     * @brief Whether the database has already computed a new paramter point
     *        requiring a residual.
     */
    bool d_next_parameter_point_computed;

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
     * @brief The rank of the given processor.
     */
    int d_rank;

    /**
     * @brief Random engine used to generate subsets.
     */
    std::default_random_engine rng;
};

// Create a greedy residual point.
struct GreedyResidualPoint createGreedyResidualPoint(Vector* point, Vector* localROM);
struct GreedyResidualPoint createGreedyResidualPoint(Vector* point, std::shared_ptr<Vector>& localROM);

// Given a a vector/double, find the nearest point in a domain.
Vector getNearestPoint(std::vector<Vector> paramPoints, Vector point);
double getNearestPoint(std::vector<double> paramPoints, double point);
int getNearestPointIndex(std::vector<Vector> paramPoints, Vector point);
int getNearestPointIndex(std::vector<double> paramPoints, double point);

}

#endif
