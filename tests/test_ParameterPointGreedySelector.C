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
#include "../ParameterPointGreedySelector.h"
#define _USE_MATH_DEFINES
#include <cmath>

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
  SUCCEED();
}

TEST(ParameterPointGreedySelectorSerialTest, Test_ParameterPointGreedySelectorCentroid)
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
    CAROM::ParameterPointGreedySelector caromGreedySelector(paramPoints, 1, 1, 3, 4, true, 1, true);

    int wrongOrder = caromGreedySelector.computePointRequiringResidual();
    EXPECT_EQ(wrongOrder, -1);

    int nextPointToSample = caromGreedySelector.computeNextSampleParameterPoint();
    EXPECT_EQ(nextPointToSample, 2);
    EXPECT_EQ(paramPoints[nextPointToSample], 3.0);

    wrongOrder = caromGreedySelector.computeNextSampleParameterPoint();
    EXPECT_EQ(wrongOrder, -1);
    EXPECT_DEATH(caromGreedySelector.setPointResidualError(100.0, d_rank, d_num_procs), ".*");

    // ERRORS: [INF, INF, 0, INF, INF]

    int firstPoint = caromGreedySelector.computePointRequiringResidual();
    EXPECT_EQ(firstPoint, 0);

    wrongOrder = caromGreedySelector.computeNextSampleParameterPoint();
    EXPECT_EQ(wrongOrder, -1);
    wrongOrder = caromGreedySelector.computePointRequiringResidual();
    EXPECT_EQ(wrongOrder, -1);

    // EXPECT_DEATH(caromGreedySelector.setPointResidualError(-1.0, d_rank, d_num_procs), ".*");
    caromGreedySelector.setPointResidualError(100.0, d_rank, d_num_procs);

    wrongOrder = caromGreedySelector.computeNextSampleParameterPoint();
    EXPECT_EQ(wrongOrder, -1);

    int secondPoint = caromGreedySelector.computePointRequiringResidual();
    EXPECT_EQ(secondPoint, 1);
    caromGreedySelector.setPointResidualError(50.0, d_rank, d_num_procs);
    int thirdPoint = caromGreedySelector.computePointRequiringResidual();
    EXPECT_EQ(thirdPoint, 3);
    caromGreedySelector.setPointResidualError(30.0, d_rank, d_num_procs);
    nextPointToSample = caromGreedySelector.computeNextSampleParameterPoint();
    EXPECT_EQ(paramPoints[nextPointToSample], paramPoints[firstPoint]);

    // ERRORS: [0, 50, 0, 30, INF, INF]

    firstPoint = caromGreedySelector.computePointRequiringResidual();
    EXPECT_EQ(firstPoint, 1);
    caromGreedySelector.setPointResidualError(0.9, d_rank, d_num_procs);
    secondPoint = caromGreedySelector.computePointRequiringResidual();
    EXPECT_EQ(secondPoint, 3);
    caromGreedySelector.setPointResidualError(25.0, d_rank, d_num_procs);
    thirdPoint = caromGreedySelector.computePointRequiringResidual();
    EXPECT_EQ(thirdPoint, 4);
    caromGreedySelector.setPointResidualError(0.9, d_rank, d_num_procs);
    nextPointToSample = caromGreedySelector.computeNextSampleParameterPoint();
    EXPECT_EQ(paramPoints[nextPointToSample], paramPoints[secondPoint]);

    // ERRORS: [0, 0.9, 0, 0.9, 0.9, INF]

    firstPoint = caromGreedySelector.computePointRequiringResidual();
    caromGreedySelector.setPointResidualError(0.9, d_rank, d_num_procs);
    EXPECT_EQ(firstPoint, 1);
    secondPoint = caromGreedySelector.computePointRequiringResidual();
    EXPECT_EQ(secondPoint, 5);
    caromGreedySelector.setPointResidualError(0.9, d_rank, d_num_procs);

    // ERRORS: [0, 0.9, 0, 0, 0.9, 0.9]

    nextPointToSample = caromGreedySelector.computeNextSampleParameterPoint();
    EXPECT_EQ(nextPointToSample, -1);
}

TEST(ParameterPointGreedySelectorSerialTest, Test_ParameterPointGreedySelector)
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
    CAROM::ParameterPointGreedySelector caromGreedySelector(paramPoints, 1, 1, 3, 4, false, 1, true);

    int nextPointToSample = caromGreedySelector.computeNextSampleParameterPoint();
    EXPECT_EQ(nextPointToSample, 2);
    EXPECT_EQ(paramPoints[nextPointToSample], 3.0);

    int closestROMIndex = caromGreedySelector.obtainNearestROMIndex(0);
    EXPECT_EQ(closestROMIndex, 2);

    // ERRORS: [INF, INF, 0, INF, INF]

    int firstPoint = caromGreedySelector.computePointRequiringResidual();
    EXPECT_EQ(firstPoint, 0);
    caromGreedySelector.setPointResidualError(100.0, d_rank, d_num_procs);

    int secondPoint = caromGreedySelector.computePointRequiringResidual();
    EXPECT_EQ(secondPoint, 1);
    caromGreedySelector.setPointResidualError(50.0, d_rank, d_num_procs);
    int thirdPoint = caromGreedySelector.computePointRequiringResidual();
    EXPECT_EQ(thirdPoint, 3);
    caromGreedySelector.setPointResidualError(30.0, d_rank, d_num_procs);
    nextPointToSample = caromGreedySelector.computeNextSampleParameterPoint();
    EXPECT_EQ(paramPoints[nextPointToSample], paramPoints[firstPoint]);

    closestROMIndex = caromGreedySelector.obtainNearestROMIndex(0);
    EXPECT_EQ(closestROMIndex, 2);

    // ERRORS: [0, 50, 0, 30, INF, INF]

    firstPoint = caromGreedySelector.computePointRequiringResidual();
    EXPECT_EQ(firstPoint, 1);
    caromGreedySelector.setPointResidualError(0.9, d_rank, d_num_procs);
    secondPoint = caromGreedySelector.computePointRequiringResidual();
    EXPECT_EQ(secondPoint, 3);
    caromGreedySelector.setPointResidualError(25.0, d_rank, d_num_procs);
    thirdPoint = caromGreedySelector.computePointRequiringResidual();
    EXPECT_EQ(thirdPoint, 4);
    caromGreedySelector.setPointResidualError(0.9, d_rank, d_num_procs);
    nextPointToSample = caromGreedySelector.computeNextSampleParameterPoint();
    EXPECT_EQ(paramPoints[nextPointToSample], paramPoints[secondPoint]);

    closestROMIndex = caromGreedySelector.obtainNearestROMIndex(1);
    EXPECT_EQ(closestROMIndex, 0);

    EXPECT_DEATH(caromGreedySelector.obtainNearestROMIndex(-1), ".*");
    EXPECT_DEATH(caromGreedySelector.obtainNearestROMIndex(10), ".*");
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  return RUN_ALL_TESTS();
}
#else // #ifndef CAROM_HAS_GTEST
int main()
{
  std::cout << "libROM was compiled without Google Test support, so unit "
	    << "tests have been disabled. To enable unit tests, compile "
	    << "libROM with Google Test support." << std::endl;
}
#endif // #endif CAROM_HAS_GTEST
