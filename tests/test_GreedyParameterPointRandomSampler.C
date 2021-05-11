/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: This source file is a test runner that uses the Google Test
// Framework to run unit tests on the CAROM::GNAT class.

#include <iostream>

#ifdef CAROM_HAS_GTEST
#include<gtest/gtest.h>
#include <mpi.h>
#include "../GreedyParameterPointRandomSampler.h"
#define _USE_MATH_DEFINES
#include <cmath>

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
    SUCCEED();
}

TEST(GreedyParameterPointRandomSamplerSerialTest, Test_GreedyParameterPointRandomSamplerCentroid)
{
    // Get the rank of this process, and the number of processors.
    int mpi_init, d_rank, d_num_procs;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

    CAROM::GreedyParameterPointRandomSampler caromGreedySampler(0.5, 2.5, 5, false, 0.1, 1, 1, 3, 4, false, "", "", true, 1, true);

    std::vector<CAROM::Vector> paramDomain = caromGreedySampler.getParameterPointDomain();

    std::shared_ptr<CAROM::Vector> nextPointToSample = caromGreedySampler.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), paramDomain[3].item(0));

    // ERRORS: [INF, INF, 0, INF, INF]
    caromGreedySampler.getNextPointRequiringRelativeError();
    caromGreedySampler.setPointRelativeError(100.0);

    struct CAROM::GreedyErrorIndicatorPoint localPoint = caromGreedySampler.getNextPointRequiringErrorIndicator();
    EXPECT_EQ(localPoint.point.get()->item(0), paramDomain[3].item(0));
    caromGreedySampler.setPointErrorIndicator(1.0, 1);

    struct CAROM::GreedyErrorIndicatorPoint firstPoint = caromGreedySampler.getNextPointRequiringErrorIndicator();
    EXPECT_EQ(firstPoint.point.get()->item(0), paramDomain[0].item(0));
    caromGreedySampler.setPointErrorIndicator(100.0, 1);
    struct CAROM::GreedyErrorIndicatorPoint secondPoint = caromGreedySampler.getNextPointRequiringErrorIndicator();
    EXPECT_EQ(secondPoint.point.get()->item(0), paramDomain[1].item(0));
    caromGreedySampler.setPointErrorIndicator(50.0, 1);
    struct CAROM::GreedyErrorIndicatorPoint thirdPoint = caromGreedySampler.getNextPointRequiringErrorIndicator();
    EXPECT_EQ(thirdPoint.point.get()->item(0), paramDomain[2].item(0));
    caromGreedySampler.setPointErrorIndicator(30.0, 1);

    nextPointToSample = caromGreedySampler.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), firstPoint.point.get()->item(0));

    // ERRORS: [0, 50, 0, 30, INF, INF]
    caromGreedySampler.getNextPointRequiringRelativeError();
    caromGreedySampler.setPointRelativeError(100.0);

    caromGreedySampler.getNextPointRequiringErrorIndicator();
    EXPECT_EQ(firstPoint.point.get()->item(0), paramDomain[0].item(0));
    caromGreedySampler.setPointErrorIndicator(35.0, 1);
    nextPointToSample = caromGreedySampler.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), paramDomain[1].item(0));

    // ERRORS: [0, 0, 0, 30, 35, INF]
    caromGreedySampler.getNextPointRequiringRelativeError();
    caromGreedySampler.setPointRelativeError(100.0);

    firstPoint = caromGreedySampler.getNextPointRequiringErrorIndicator();
    caromGreedySampler.setPointErrorIndicator(0.3, 1);
    EXPECT_EQ(firstPoint.point.get()->item(0), paramDomain[2].item(0));
    secondPoint = caromGreedySampler.getNextPointRequiringErrorIndicator();
    caromGreedySampler.setPointErrorIndicator(0.4, 1);
    EXPECT_EQ(secondPoint.point.get()->item(0), paramDomain[4].item(0));
    nextPointToSample = caromGreedySampler.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), paramDomain[4].item(0));

    // ERRORS: [0, 0, 0, 30, 0, 0.3]
    caromGreedySampler.getNextPointRequiringRelativeError();
    caromGreedySampler.setPointRelativeError(100.0);

    firstPoint = caromGreedySampler.getNextPointRequiringErrorIndicator();
    EXPECT_EQ(firstPoint.point.get(), nullptr);
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
#else // #ifndef CAROM_HAS_GTEST
int main()
{
    std::cout << "libROM was compiled without Google Test support, so unit "
              << "tests have been disabled. To enable unit tests, compile "
              << "libROM with Google Test support." << std::endl;
}
#endif // #endif CAROM_HAS_GTEST
