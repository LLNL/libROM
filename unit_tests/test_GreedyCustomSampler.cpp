/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
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
#include "algo/greedy/GreedyCustomSampler.h"
#define _USE_MATH_DEFINES
#include <cmath>

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
    SUCCEED();
}

TEST(GreedyCustomSamplerSerialTest, Test_GreedyCustomSamplerCentroid)
{
    // Get the rank of this process, and the number of processors.
    int mpi_init, d_rank, d_num_procs;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

    std::vector<double> paramPoints = {1.0, 3.0, 6.0};
    CAROM::GreedyCustomSampler caromGreedySampler(paramPoints, false, 0.1, 1, 1, 2,
            3, "", "", true, 1, true);

    std::shared_ptr<CAROM::Vector> nextPointToSample =
        caromGreedySampler.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), 3.0);

    caromGreedySampler.getNextPointRequiringRelativeError();
    caromGreedySampler.setPointRelativeError(100.0);

    struct CAROM::GreedyErrorIndicatorPoint localPoint =
        caromGreedySampler.getNextPointRequiringErrorIndicator();
    EXPECT_EQ(localPoint.point.get()->item(0), 3.0);
    caromGreedySampler.setPointErrorIndicator(1.0, 1);

    struct CAROM::GreedyErrorIndicatorPoint firstPoint =
        caromGreedySampler.getNextPointRequiringErrorIndicator();
    EXPECT_EQ(firstPoint.point.get()->item(0), 1.0);
    caromGreedySampler.setPointErrorIndicator(100.0, 1);
    struct CAROM::GreedyErrorIndicatorPoint secondPoint =
        caromGreedySampler.getNextPointRequiringErrorIndicator();
    EXPECT_EQ(secondPoint.point.get()->item(0), 6.0);
    caromGreedySampler.setPointErrorIndicator(50.0, 1);

    nextPointToSample = caromGreedySampler.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), firstPoint.point.get()->item(0));

    caromGreedySampler.getNextPointRequiringRelativeError();
    caromGreedySampler.setPointRelativeError(100.0);

    firstPoint = caromGreedySampler.getNextPointRequiringErrorIndicator();
    nextPointToSample = caromGreedySampler.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), paramPoints[2]);

    caromGreedySampler.getNextPointRequiringRelativeError();
    caromGreedySampler.setPointRelativeError(100.0);

    firstPoint = caromGreedySampler.getNextPointRequiringErrorIndicator();
}

TEST(GreedyCustomSamplerSerialTest,
     Test_GreedyCustomSamplerMultiDimensionalCentroid)
{
    // Get the rank of this process, and the number of processors.
    int mpi_init, d_rank, d_num_procs;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

    CAROM::Vector item0(2, false);
    item0.item(0) = 1.0;
    item0.item(1) = 11.0;

    CAROM::Vector item1(2, false);
    item1.item(0) = 3.0;
    item1.item(1) = 13.0;

    CAROM::Vector item2(2, false);
    item2.item(0) = 6.0;
    item2.item(1) = 16.0;

    std::vector<CAROM::Vector> paramPoints = {item0, item1, item2};
    CAROM::GreedyCustomSampler caromGreedySampler(paramPoints, false, 0.1, 1, 1, 2,
            3, "", "", true, 1, true);

    std::shared_ptr<CAROM::Vector> nextPointToSample =
        caromGreedySampler.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), 3.0);
    EXPECT_EQ(nextPointToSample.get()->item(1), 13.0);

    caromGreedySampler.getNextPointRequiringRelativeError();
    caromGreedySampler.setPointRelativeError(100.0);

    struct CAROM::GreedyErrorIndicatorPoint localPoint =
        caromGreedySampler.getNextPointRequiringErrorIndicator();
    EXPECT_EQ(localPoint.point.get()->item(0), 3.0);
    EXPECT_EQ(localPoint.point.get()->item(1), 13.0);
    caromGreedySampler.setPointErrorIndicator(1.0, 1);

    struct CAROM::GreedyErrorIndicatorPoint firstPoint =
        caromGreedySampler.getNextPointRequiringErrorIndicator();
    EXPECT_EQ(firstPoint.point.get()->item(0), 1.0);
    EXPECT_EQ(firstPoint.point.get()->item(1), 11.0);
    caromGreedySampler.setPointErrorIndicator(100.0, 1);
    struct CAROM::GreedyErrorIndicatorPoint secondPoint =
        caromGreedySampler.getNextPointRequiringErrorIndicator();
    EXPECT_EQ(secondPoint.point.get()->item(0), 6.0);
    EXPECT_EQ(secondPoint.point.get()->item(1), 16.0);
    caromGreedySampler.setPointErrorIndicator(50.0, 1);

    nextPointToSample = caromGreedySampler.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), firstPoint.point.get()->item(0));
    EXPECT_EQ(nextPointToSample.get()->item(1), firstPoint.point.get()->item(1));

    caromGreedySampler.getNextPointRequiringRelativeError();
    caromGreedySampler.setPointRelativeError(100.0);

    firstPoint = caromGreedySampler.getNextPointRequiringErrorIndicator();
    nextPointToSample = caromGreedySampler.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), item2.item(0));
    EXPECT_EQ(nextPointToSample.get()->item(1), item2.item(1));

    caromGreedySampler.getNextPointRequiringRelativeError();
    caromGreedySampler.setPointRelativeError(100.0);

    firstPoint = caromGreedySampler.getNextPointRequiringErrorIndicator();
}

TEST(GreedyCustomSamplerSerialTest, Test_GreedySaveAndLoad)
{
    // Get the rank of this process, and the number of processors.
    int mpi_init, d_rank, d_num_procs;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

    std::vector<double> paramPoints = {1.0, 2.0, 3.0, 99.0, 100.0, 101.0};
    CAROM::GreedyCustomSampler caromGreedySampler(paramPoints, false, 0.1, 1, 1, 3,
            4, "", "", false, 1, true);
    caromGreedySampler.save("greedy_test");

    CAROM::GreedyCustomSampler caromGreedySamplerLoad("greedy_test");
    caromGreedySamplerLoad.save("greedy_test_LOAD");

    CAROM::Vector pointToFindNearestROM(1, false);
    pointToFindNearestROM.item(0) = 1.0;

    std::shared_ptr<CAROM::Vector> closestROM = caromGreedySampler.getNearestROM(
                pointToFindNearestROM);
    std::shared_ptr<CAROM::Vector> closestROMLoad =
        caromGreedySamplerLoad.getNearestROM(pointToFindNearestROM);
    EXPECT_EQ(closestROM.get(), closestROMLoad.get());

    std::shared_ptr<CAROM::Vector> nextPointToSample =
        caromGreedySampler.getNextParameterPoint();
    std::shared_ptr<CAROM::Vector> nextPointToSampleLoad =
        caromGreedySamplerLoad.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0),
              nextPointToSampleLoad.get()->item(0));
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
