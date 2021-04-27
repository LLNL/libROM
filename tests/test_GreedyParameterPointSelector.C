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
#include "../GreedyParameterPointSelector.h"
#define _USE_MATH_DEFINES
#include <cmath>

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
    SUCCEED();
}

TEST(GreedyParameterPointSelectorSerialTest, Test_GreedyParameterPointSelectorCentroid)
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
    CAROM::GreedyParameterPointSelector caromGreedySelector(paramPoints, false, 1, 1, 3, 4, "", "", true, 1, true);

    std::shared_ptr<CAROM::Vector> nextPointToSample = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), 3.0);

    // ERRORS: [INF, INF, 0, INF, INF]

    struct CAROM::GreedyResidualPoint firstPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(firstPoint.point.get()->item(0), 1.0);

    caromGreedySelector.setPointResidual(100.0, 1);

    struct CAROM::GreedyResidualPoint secondPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(secondPoint.point.get()->item(0), 2.0);
    caromGreedySelector.setPointResidual(50.0, 1);
    struct CAROM::GreedyResidualPoint thirdPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(thirdPoint.point.get()->item(0), 99.0);
    caromGreedySelector.setPointResidual(30.0, 1);
    nextPointToSample = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), firstPoint.point.get()->item(0));

    // ERRORS: [0, 50, 0, 30, INF, INF]

    firstPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(firstPoint.point.get()->item(0), 100.0);
    caromGreedySelector.setPointResidual(35.0, 1);
    nextPointToSample = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), paramPoints[1]);

    // ERRORS: [0, 0, 0, 30, 35, INF]

    firstPoint = caromGreedySelector.getNextPointRequiringResidual();
    caromGreedySelector.setPointResidual(0.3, 1);
    EXPECT_EQ(firstPoint.point.get()->item(0), 101.0);
    nextPointToSample = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), paramPoints[4]);

    // ERRORS: [0, 0, 0, 30, 0, 0.3]

    firstPoint = caromGreedySelector.getNextPointRequiringResidual();
    caromGreedySelector.setPointResidual(0.3, 1);
    EXPECT_EQ(firstPoint.point.get()->item(0), 99.0);
}

TEST(GreedyParameterPointSelectorSerialTest, Test_GreedyParameterPointSelector)
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
    CAROM::GreedyParameterPointSelector caromGreedySelector(paramPoints, false, 1, 1, 3, 4, "", "", false, 1, true);

    CAROM::Vector pointToFindNearestROM(1, false);
    pointToFindNearestROM.item(0) = 1.0;

    std::shared_ptr<CAROM::Vector> closestROM = caromGreedySelector.getNearestROM(pointToFindNearestROM);
    EXPECT_EQ(closestROM.get(), nullptr);

    std::shared_ptr<CAROM::Vector> nextPointToSample = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), 3.0);

    closestROM = caromGreedySelector.getNearestROM(pointToFindNearestROM);
    EXPECT_EQ(closestROM.get()->item(0), 3.0);

    // ERRORS: [INF, INF, 0, INF, INF]

    struct CAROM::GreedyResidualPoint firstPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(firstPoint.point.get()->item(0), 1.0);
    caromGreedySelector.setPointResidual(100.0, 1);
    struct CAROM::GreedyResidualPoint secondPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(secondPoint.point.get()->item(0), 2.0);
    caromGreedySelector.setPointResidual(50.0, 1);
    struct CAROM::GreedyResidualPoint thirdPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(thirdPoint.point.get()->item(0), 99.0);
    caromGreedySelector.setPointResidual(30.0, 1);

    closestROM = caromGreedySelector.getNearestROM(pointToFindNearestROM);
    EXPECT_EQ(closestROM.get()->item(0), 3.0);
    nextPointToSample = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), firstPoint.point.get()->item(0));

    // ERRORS: [0, 50, 0, 30, INF, INF]

    firstPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(firstPoint.point.get()->item(0), 100.0);
    caromGreedySelector.setPointResidual(0.9, 1);
    nextPointToSample = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), paramPoints[1]);

    pointToFindNearestROM.item(0) = 2.0;

    closestROM = caromGreedySelector.getNearestROM(pointToFindNearestROM);
    EXPECT_EQ(closestROM.get()->item(0), 2.0);

    pointToFindNearestROM.item(0) = 100.0;

    closestROM = caromGreedySelector.getNearestROM(pointToFindNearestROM);
    EXPECT_EQ(closestROM.get()->item(0), 3.0);

}

TEST(GreedyParameterPointSelectorSerialTest, Test_GreedyParameterPointSaveAndLoad)
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
    CAROM::GreedyParameterPointSelector caromGreedySelector(paramPoints, false, 1, 1, 3, 4, "", "", false, 1, true);

    caromGreedySelector.save("greedy_test");

    CAROM::GreedyParameterPointSelector caromGreedySelectorLoad("greedy_test");
    caromGreedySelectorLoad.save("greedy_test_LOAD");

    CAROM::Vector pointToFindNearestROM(1, false);
    pointToFindNearestROM.item(0) = 1.0;

    std::shared_ptr<CAROM::Vector> closestROM = caromGreedySelector.getNearestROM(pointToFindNearestROM);
    std::shared_ptr<CAROM::Vector> closestROMLoad = caromGreedySelectorLoad.getNearestROM(pointToFindNearestROM);
    EXPECT_EQ(closestROM.get(), closestROMLoad.get());

    std::shared_ptr<CAROM::Vector> nextPointToSample = caromGreedySelector.getNextParameterPoint();
    std::shared_ptr<CAROM::Vector> nextPointToSampleLoad = caromGreedySelectorLoad.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), nextPointToSampleLoad.get()->item(0));

    closestROM = caromGreedySelector.getNearestROM(pointToFindNearestROM);
    closestROMLoad = caromGreedySelectorLoad.getNearestROM(pointToFindNearestROM);
    EXPECT_EQ(closestROM.get()->item(0), closestROMLoad.get()->item(0));

    // ERRORS: [INF, INF, 0, INF, INF]

    struct CAROM::GreedyResidualPoint firstPoint = caromGreedySelector.getNextPointRequiringResidual();
    struct CAROM::GreedyResidualPoint firstPointLoad = caromGreedySelectorLoad.getNextPointRequiringResidual();
    EXPECT_EQ(firstPoint.point.get()->item(0), firstPointLoad.point.get()->item(0));
    caromGreedySelector.setPointResidual(100.0, 1);
    caromGreedySelectorLoad.setPointResidual(100.0, 1);
    struct CAROM::GreedyResidualPoint secondPoint = caromGreedySelector.getNextPointRequiringResidual();
    struct CAROM::GreedyResidualPoint secondPointLoad = caromGreedySelectorLoad.getNextPointRequiringResidual();
    EXPECT_EQ(secondPoint.point.get()->item(0), secondPointLoad.point.get()->item(0));
    caromGreedySelector.setPointResidual(50.0, 1);
    caromGreedySelectorLoad.setPointResidual(50.0, 1);
    struct CAROM::GreedyResidualPoint thirdPoint = caromGreedySelector.getNextPointRequiringResidual();
    struct CAROM::GreedyResidualPoint thirdPointLoad = caromGreedySelectorLoad.getNextPointRequiringResidual();
    EXPECT_EQ(thirdPoint.point.get()->item(0), thirdPointLoad.point.get()->item(0));
    caromGreedySelector.setPointResidual(30.0, 1);
    caromGreedySelectorLoad.setPointResidual(30.0, 1);

    closestROM = caromGreedySelector.getNearestROM(pointToFindNearestROM);
    closestROMLoad = caromGreedySelectorLoad.getNearestROM(pointToFindNearestROM);
    EXPECT_EQ(closestROM.get()->item(0), closestROMLoad.get()->item(0));
    nextPointToSample = caromGreedySelector.getNextParameterPoint();
    nextPointToSampleLoad = caromGreedySelectorLoad.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample.get()->item(0), nextPointToSampleLoad.get()->item(0));

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
