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
class IOptimizable
{
public:
    struct Constraints
    {
        Constraints(double lower = 0.0, double upper = 1.0,
                    bool isConstrained = false) :
            lower(lower),
            upper(upper),
            isConstrained(isConstrained)
        {

        }

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

        double lower;
        double upper;
        bool isConstrained;
    };

    virtual double EvaluateCost(std::vector<double> inputs) const = 0;
    virtual unsigned int NumberOfParameters() const = 0;
    virtual std::vector<Constraints> GetConstraints() const = 0;
    virtual ~IOptimizable() {}
};

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
     * @param[in] shouldCheckConstraints Should constraints bee checked on for each new candidate.
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
    void Optimize(int min_iterations, int max_iterations, double cost_tolerance,
                  bool verbose = true);

private:
    bool CheckConstraints(std::vector<double> agent);

    void InitPopulation();

    void SelectionAndCrossing();

    std::vector<double> GetBestAgent() const;

    double GetBestCost() const;

    std::vector<std::pair<std::vector<double>, double>> GetPopulationWithCosts()
            const;

    void PrintPopulation() const;

    const IOptimizable& m_cost;
    unsigned int m_populationSize;
    double m_F;
    double m_CR;

    unsigned int m_numberOfParameters;

    bool m_shouldCheckConstraints;

    std::function<void(const DifferentialEvolution&)> m_callback;
    std::function<bool(const DifferentialEvolution&)> m_terminationCondition;

    std::default_random_engine m_generator;
    std::vector<std::vector<double>> m_population;

    std::vector<double> m_minCostPerAgent;

    std::vector<IOptimizable::Constraints> m_constraints;

    int m_bestAgentIndex;
    double m_minCost;

    static constexpr double g_defaultLowerConstraint =
        -std::numeric_limits<double>::infinity();
    static constexpr double g_defaultUpperConstarint =
        std::numeric_limits<double>::infinity();
};
}

#endif
