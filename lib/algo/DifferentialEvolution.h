/******************************************************************************
 *
 * Author: Milos Stojanovic Stojke (milsto)
 * Adapted by Kevin Huynh from DifferentialEvolution.h from milsto's Github repo,
 * https://github.com/milsto/differential-evolution, permission granted by
 * milsto/differential-evolution's MIT license.
 *
 * SPDX-License-Identifier: (MIT)
 *
 *****************************************************************************/

// Description: Implementation for the differential evolution algorithm.
//              The implemented differential evolution algorithm is derived from
//              Storn, R., Price, K. Differential Evolution – A Simple and
//              Efficient Heuristic for global Optimization over Continuous
//              Spaces. Journal of Global Optimization 11, 341–359 (1997).
//              https://doi.org/10.1023/A:1008202821328

#ifndef included_DifferentialEvolution_h
#define included_DifferentialEvolution_h

#include <vector>
#include <functional>
#include <random>
#include <utility>
#include <memory>
#include <limits>

namespace CAROM
{
/**
 * Class IOptimizable is a simple interface for defining the constrained optimization problem.
 */
class IOptimizable
{
public:
    /**
     * @brief Constraints to be fulfilled by the variable candidates.
     */
    struct Constraints
    {
        /**
         * @brief Constructor.
         *
         * @param[in] lower           Lower bound for constraint
         * @param[in] upper           Upper bound for constraint
         * @param[in] isConstrained   Whether to check the constraints
         */
        Constraints(double lower = 0.0, double upper = 1.0,
                    bool isConstrained = false) :
            lower(lower),
            upper(upper),
            isConstrained(isConstrained)
        {

        }

        /**
         * @brief Check whether the current variable candidate fulfills the constraints.
         */
        bool Check(double candidate)
        {
            if (isConstrained)
            {
                if (candidate <= upper && candidate >= lower)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return true;
            }
        }

        /**
         * @brief Lower bound for constraint
         */
        double lower;

        /**
         * @brief Upper bound for constraint
         */
        double upper;

        /**
         * @brief Whether to check the constraints
         */
        bool isConstrained;
    };

    /**
     * @brief Evaluate the cost function with the current set of inputs.
     */
    virtual double EvaluateCost(std::vector<double> inputs) const = 0;

    /**
     * @brief Return the number of parameters.
     */
    virtual unsigned int NumberOfParameters() const = 0;

    /**
     * @brief Return the list of constraints.
     */
    virtual std::vector<Constraints> GetConstraints() const = 0;

    /**
     * @brief Destructor.
     */
    virtual ~IOptimizable() {}
};

/**
 * Class DifferentialEvolution is a general purpose black box optimizer for the class IOptimizable.
 */
class DifferentialEvolution
{
public:
    /**
     * @brief Constructor.
     *
     * @param[in] costFunction           Cost function to minimize
     * @param[in] populationSize         Number of agents in each optimization step
     * @param[in] F                      The differential weight.
     * @param[in] CR                     The crossover probability.
     * @param[in] randomSeed             Set random seed to a fix value to have
     *                                   repeatable (non stochastic) experiments
     * @param[in] shouldCheckConstraints Should constraints be checked on for each new candidate.
     *                                   This check may be turned off to
     *                                   increase performance if the cost function is defined
     *                                   and has no local minimum outside of the constraints.
     * @param[in] callback               Optional callback to be called after each optimization
     *                                   iteration has finished.
     * @param[in] terminationCondition   Optional termination condition callback to be called
     *                                   after each optimization iteration has finished.
     */
    DifferentialEvolution(const IOptimizable& costFunction,
                          unsigned int populationSize,
                          double F = 0.8,
                          double CR = 0.9,
                          int randomSeed = 1,
                          bool shouldCheckConstraints = true,
                          std::function<void(const DifferentialEvolution&)> callback = nullptr,
                          std::function<bool(const DifferentialEvolution&)> terminationCondition =
                              nullptr);

    /**
     * @brief Constructor.
     *
     * @param[in] min_iterations         The minimum number of iterations to run.
     *                                   Cost tolerance is not checked until then.
     * @param[in] max_iterations         The maximum number of iterations to run.
     * @param[in] cost_tolerance         The cost tolerance to determine convergence.
     * @param[in] verbose                Verbosity.
     */
    std::vector<double> Optimize(int min_iterations, int max_iterations,
                                 double cost_tolerance,
                                 bool verbose = true);

private:

    /**
     * @brief Check that the agent meets the constraints.
     */
    bool CheckConstraints(std::vector<double> agent);

    /**
     * @brief Initialize the population at the first iteration.
     */
    void InitPopulation();

    /**
     * @brief Select and cross agents during each iteration.
     */
    void SelectionAndCrossing();

    /**
     * @brief Get current best agent.
     */
    std::vector<double> GetBestAgent() const;

    /**
     * @brief Get current best cost.
     */
    double GetBestCost() const;

    /**
     * @brief Get population and its associated cost.
     */
    std::vector<std::pair<std::vector<double>, double>> GetPopulationWithCosts()
            const;

    /**
     * @brief Print population information.
     */
    void PrintPopulation() const;

    /**
     * @brief Cost function
     */
    const IOptimizable& m_cost;

    /**
     * @brief Population size
     */
    unsigned int m_populationSize;

    /**
     * @brief Differential weight
     */
    double m_F;

    /**
     * @brief Crossover probability
     */
    double m_CR;

    /**
     * @brief Number of parameters
     */
    unsigned int m_numberOfParameters;

    /**
     * @brief Whether to check constraints.
     */
    bool m_shouldCheckConstraints;

    /**
     * @brief Optional callback function.
     */
    std::function<void(const DifferentialEvolution&)> m_callback;

    /**
     * @brief Optional termination condition callback function.
     */
    std::function<bool(const DifferentialEvolution&)> m_terminationCondition;

    /**
     * @brief Random number generator
     */
    std::default_random_engine m_generator;

    /**
     * @brief Container holding populations
     */
    std::vector<std::vector<double>> m_population;

    /**
     * @brief Container holding minimum cost per agent
     */
    std::vector<double> m_minCostPerAgent;

    /**
     * @brief Container holding constraints
     */
    std::vector<IOptimizable::Constraints> m_constraints;

    /**
     * @brief Index of the current best agent
     */
    int m_bestAgentIndex;

    /**
     * @brief Current minimum cost
     */
    double m_minCost;

    /**
     * @brief Lower and upper constraint limits
     */
    static constexpr double g_defaultLowerConstraint =
        -std::numeric_limits<double>::infinity();
    static constexpr double g_defaultUpperConstraint =
        std::numeric_limits<double>::infinity();

    /**
     * @brief The rank of the process this object belongs to.
     */
    int d_rank;
};
}

#endif
