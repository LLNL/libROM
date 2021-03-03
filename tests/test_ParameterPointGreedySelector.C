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

    std::vector<double> paramPoints = {1.0, 2.0, 3.0, 20.0};
    CAROM::ParameterPointGreedySelector caromGreedySelector(paramPoints, 1, 1, 2, 3, true, 1);

    int nextPointToSample = caromGreedySelector.computeNextSampleParameterPoint();
    EXPECT_EQ(paramPoints[nextPointToSample], 3.0);
}

TEST(ParameterPointGreedySelectorSerialTest, Test_ParameterPointGreedySelector)
{
    std::vector<double> paramPoints = {1.0, 2.0, 3.0, 20.0};
    CAROM::ParameterPointGreedySelector caromGreedySelector(paramPoints, 1, 1, 2, 3, true, 1);

    EXPECT_EQ(paramPoints[nextPointToSample], 2.0);
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
