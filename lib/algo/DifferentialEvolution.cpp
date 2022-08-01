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

// Description: Implementation of the differential evolution algorithm.

#include "DifferentialEvolution.h"
#include "utils/Utilities.h"
#include "mpi.h"

#include <iostream>

namespace CAROM {

DifferentialEvolution::DifferentialEvolution(const IOptimizable& costFunction,
        unsigned int populationSize,
        double F,
        double CR,
        int randomSeed,
        bool shouldCheckConstraints,
        std::function<void(const DifferentialEvolution&)> callback,
        std::function<bool(const DifferentialEvolution&)> terminationCondition) :
    m_cost(costFunction),
    m_populationSize(populationSize),
    m_F(F),
    m_CR(CR),
    m_bestAgentIndex(0),
    m_minCost(-std::numeric_limits<double>::infinity()),
    m_shouldCheckConstraints(shouldCheckConstraints),
    m_callback(callback),
    m_terminationCondition(terminationCondition)
{
    // Get the rank of this process, and the number of processors.
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);

    m_generator.seed(randomSeed);
    CAROM_VERIFY(m_populationSize >= 4);

    m_numberOfParameters = m_cost.NumberOfParameters();

    m_population.resize(populationSize);
    for (auto& agent : m_population)
    {
        agent.resize(m_numberOfParameters);
    }

    m_minCostPerAgent.resize(m_populationSize);

    m_constraints = costFunction.GetConstraints();
}

std::vector<double>
DifferentialEvolution::Optimize(int min_iterations, int max_iterations,
                                double cost_tolerance, bool verbose)
{
    CAROM_VERIFY(min_iterations <= max_iterations);
    CAROM_VERIFY(min_iterations >= 0 && max_iterations > 0);

    std::vector<double> optimal_parameters(m_numberOfParameters);

    InitPopulation();

    double prevMinCost;

    // Optimization loop
    for (int i = 0; i < max_iterations; i++)
    {
        // Optimization step
        prevMinCost = m_minCost;
        SelectionAndCrossing();

        if (d_rank == 0 && verbose)
        {
            std::cout << "Iteration: " << i << "\t\t";
            std::cout << "Current minimal cost: " << m_minCost << "\t\t";
            std::cout << "Best agent: ";
            for (int j = 0; j < m_numberOfParameters; j++)
            {
                std::cout << m_population[m_bestAgentIndex][j] << " ";
            }
            std::cout << std::endl;
        }

        for (int j = 0; j < m_numberOfParameters; j++)
        {
            optimal_parameters[j] = m_population[m_bestAgentIndex][j];
        }

        if (i >= min_iterations && prevMinCost - m_minCost <= cost_tolerance)
        {
            if (d_rank == 0 && verbose)
            {
                std::cout <<
                          "Terminated due to cost tolerance condition being met." <<
                          std::endl;
            }
            return optimal_parameters;
        }

        if (m_callback)
        {
            m_callback(*this);
        }

        if (m_terminationCondition)
        {
            if (m_terminationCondition(*this))
            {
                if (d_rank == 0 && verbose)
                {
                    std::cout <<
                              "Terminated due to positive evaluation of the termination condition." <<
                              std::endl;
                }
                return optimal_parameters;
            }
        }
    }

    if (d_rank == 0 && verbose)
    {
        std::cout << "Terminated due to exceeding total number of generations." <<
                  std::endl;
    }

    return optimal_parameters;
}

bool
DifferentialEvolution::CheckConstraints(std::vector<double> agent)
{
    for (int i = 0; i < agent.size(); i++)
    {
        if (!m_constraints[i].Check(agent[i]))
        {
            return false;
        }
    }

    return true;
}

void
DifferentialEvolution::InitPopulation()
{
    // Init population based on random sampling of the cost function
    std::shared_ptr<std::uniform_real_distribution<double>> distribution;

    for (auto& agent : m_population)
    {
        for (int i = 0; i < m_numberOfParameters; i++)
        {
            if (m_constraints[i].isConstrained)
            {
                distribution = std::make_shared<std::uniform_real_distribution<double>>
                               (std::uniform_real_distribution<double>(m_constraints[i].lower,
                                       m_constraints[i].upper));
            }
            else
            {
                distribution = std::make_shared<std::uniform_real_distribution<double>>
                               (std::uniform_real_distribution<double>(g_defaultLowerConstraint,
                                       g_defaultUpperConstraint));
            }

            agent[i] = (*distribution)(m_generator);
        }
    }

    // Initialize minimum cost, best agent and best agent index
    for (int i = 0; i < m_populationSize; i++)
    {
        m_minCostPerAgent[i] = m_cost.EvaluateCost(m_population[i]);

        if (m_minCostPerAgent[i] < m_minCost)
        {
            m_minCost = m_minCostPerAgent[i];
            m_bestAgentIndex = i;
        }
    }
}

void
DifferentialEvolution::SelectionAndCrossing()
{
    std::uniform_real_distribution<double> distribution(0, m_populationSize);

    double minCost = m_minCostPerAgent[0];
    int bestAgentIndex = 0;

    for (int x = 0; x < m_populationSize; x++)
    {
        // For x in population select 3 random agents (a, b, c) different from x
        int a = x;
        int b = x;
        int c = x;

        // Agents must be different from each other and from x
        while (a == x || b == x || c == x || a == b || a == c || b == c)
        {
            a = distribution(m_generator);
            b = distribution(m_generator);
            c = distribution(m_generator);
        }

        // Form intermediate solution z
        std::vector<double> z(m_numberOfParameters);
        for (int i = 0; i < m_numberOfParameters; i ++)
        {
            z[i] = m_population[a][i] + m_F * (m_population[b][i] - m_population[c][i]);
        }

        // Choose random R
        std::uniform_real_distribution<double> distributionParam(0,
                m_numberOfParameters);
        int R = distributionParam(m_generator);

        // Choose random r for each dimension
        std::vector<double> r(m_numberOfParameters);
        std::uniform_real_distribution<double> distributionPerX(0, 1);
        for (auto& var : r)
        {
            var = distributionPerX(m_generator);
        }

        std::vector<double> newX(m_numberOfParameters);

        // Execute crossing
        for (int i = 0; i < m_numberOfParameters; i++)
        {
            if (r[i] < m_CR || i == R)
            {
                newX[i] = z[i];
            }
            else
            {
                newX[i] = m_population[x][i];
            }
        }

        // Check if newX candidate satisfies constraints and skip it if not.
        // If agent is skipped loop iteration x is decreased so that it is ensured
        // that the population has constant size (equal to m_populationSize).
        if (m_shouldCheckConstraints && !CheckConstraints(newX))
        {
            x--;
            continue;
        }

        // Calculate new cost and decide should the newX be kept.
        double newCost = m_cost.EvaluateCost(newX);
        if (newCost < m_minCostPerAgent[x])
        {
            m_population[x] = newX;
            m_minCostPerAgent[x] = newCost;
        }

        // Track the global best agent.
        if (m_minCostPerAgent[x] < minCost)
        {
            minCost = m_minCostPerAgent[x];
            bestAgentIndex = x;
        }
    }

    m_minCost = minCost;
    m_bestAgentIndex = bestAgentIndex;
}

std::vector<double>
DifferentialEvolution::GetBestAgent() const
{
    return m_population[m_bestAgentIndex];
}

double
DifferentialEvolution::GetBestCost() const
{
    return m_minCostPerAgent[m_bestAgentIndex];
}

std::vector<std::pair<std::vector<double>, double>>
        DifferentialEvolution::GetPopulationWithCosts() const
{
    std::vector<std::pair<std::vector<double>, double>> toRet;
    for (int i = 0; i < m_populationSize; i++)
    {
        toRet.push_back(std::make_pair(m_population[i], m_minCostPerAgent[i]));
    }

    return toRet;
}

void
DifferentialEvolution::PrintPopulation() const
{
    if (d_rank == 0)
    {
        for (auto agent : m_population)
        {
            for (auto& var : agent)
            {
                std::cout << var << " ";
            }
            std::cout << std::endl;
        }
    }
}

}
