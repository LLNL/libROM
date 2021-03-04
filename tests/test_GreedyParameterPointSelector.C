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
    CAROM::GreedyParameterPointSelector caromGreedySelector(paramPoints, 1, 1, 3, 4, true, 1, true);

    int wrongOrder = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(wrongOrder, -1);

    int nextPointToSample = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample, 2);
    EXPECT_EQ(paramPoints[nextPointToSample], 3.0);

    wrongOrder = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(wrongOrder, -1);

    // ERRORS: [INF, INF, 0, INF, INF]

    int firstPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(firstPoint, 0);

    wrongOrder = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(wrongOrder, -1);
    wrongOrder = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(wrongOrder, -1);

    caromGreedySelector.setPointResidual(100.0, d_rank, d_num_procs);

    wrongOrder = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(wrongOrder, -1);

    int secondPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(secondPoint, 1);
    caromGreedySelector.setPointResidual(50.0, d_rank, d_num_procs);
    int thirdPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(thirdPoint, 3);
    caromGreedySelector.setPointResidual(30.0, d_rank, d_num_procs);
    nextPointToSample = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(paramPoints[nextPointToSample], paramPoints[firstPoint]);

    // ERRORS: [0, 50, 0, 30, INF, INF]

    firstPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(firstPoint, 1);
    caromGreedySelector.setPointResidual(0.9, d_rank, d_num_procs);
    secondPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(secondPoint, 3);
    caromGreedySelector.setPointResidual(25.0, d_rank, d_num_procs);
    thirdPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(thirdPoint, 4);
    caromGreedySelector.setPointResidual(0.9, d_rank, d_num_procs);
    nextPointToSample = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(paramPoints[nextPointToSample], paramPoints[secondPoint]);

    // ERRORS: [0, 0.9, 0, 0.9, 0.9, INF]

    firstPoint = caromGreedySelector.getNextPointRequiringResidual();
    caromGreedySelector.setPointResidual(0.9, d_rank, d_num_procs);
    EXPECT_EQ(firstPoint, 1);
    secondPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(secondPoint, 5);
    caromGreedySelector.setPointResidual(0.9, d_rank, d_num_procs);

    // ERRORS: [0, 0.9, 0, 0, 0.9, 0.9]

    nextPointToSample = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample, -1);
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
    CAROM::GreedyParameterPointSelector caromGreedySelector(paramPoints, 1, 1, 3, 4, false, 1, true);

    int closestROMIndex = caromGreedySelector.getNearestROM(0);
    EXPECT_EQ(closestROMIndex, -1);

    int nextPointToSample = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(nextPointToSample, 2);
    EXPECT_EQ(paramPoints[nextPointToSample], 3.0);

    closestROMIndex = caromGreedySelector.getNearestROM(0);
    EXPECT_EQ(closestROMIndex, 2);

    // ERRORS: [INF, INF, 0, INF, INF]

    int firstPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(firstPoint, 0);
    caromGreedySelector.setPointResidual(100.0, d_rank, d_num_procs);
    int secondPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(secondPoint, 1);
    caromGreedySelector.setPointResidual(50.0, d_rank, d_num_procs);
    int thirdPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(thirdPoint, 3);
    caromGreedySelector.setPointResidual(30.0, d_rank, d_num_procs);

    closestROMIndex = caromGreedySelector.getNearestROM(0);
    EXPECT_EQ(closestROMIndex, 2);
    nextPointToSample = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(paramPoints[nextPointToSample], paramPoints[firstPoint]);

    // ERRORS: [0, 50, 0, 30, INF, INF]

    firstPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(firstPoint, 1);
    caromGreedySelector.setPointResidual(0.9, d_rank, d_num_procs);
    secondPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(secondPoint, 3);
    caromGreedySelector.setPointResidual(0.9, d_rank, d_num_procs);
    thirdPoint = caromGreedySelector.getNextPointRequiringResidual();
    EXPECT_EQ(thirdPoint, 4);
    caromGreedySelector.setPointResidual(25.0, d_rank, d_num_procs);
    nextPointToSample = caromGreedySelector.getNextParameterPoint();
    EXPECT_EQ(paramPoints[nextPointToSample], paramPoints[thirdPoint]);

    closestROMIndex = caromGreedySelector.getNearestROM(1);
    EXPECT_EQ(closestROMIndex, 0);

    closestROMIndex = caromGreedySelector.getNearestROM(4);
    EXPECT_EQ(closestROMIndex, 2);

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
