/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: This source file is a test runner that uses the Google Test
// Framework to run unit tests on the CAROM::DEIM class.

#include <iostream>

#ifdef CAROM_HAS_GTEST
#include<gtest/gtest.h>
#include <mpi.h>
#include "../DEIM.h"
#define _USE_MATH_DEFINES
#include <cmath>

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
  SUCCEED();
}

TEST(DEIMSerialTest, Test_DEIM)
{
  double* orthonormal_mat = {-0.1735,   -0.2978,   0.4339,   0.8304,   0.0575
   -0.5919,   -0.1608,   -0.7489,    0.2185,   -0.1230,
    0.7291,    0.1574,   -0.4449,    0.4549,   -0.1970,
   -0.2783,    0.7273,    0.1992,    0.1386,   -0.5785,
   -0.1029,    0.5760,   -0.1147,    0.1910,    0.7798}
  CAROM::Matrix u(orthonormal_mat, 5, 5, false);
  double* DEIM_res = NULL;
  int* f_sampled_row;
  int* f_sampled_rows_per_proc;
  CAROM::Matrix f_basis_sampled_inv;
  DEIM(u, 5, f_sampled_row, f_sampled_rows_per_proc, f_basis_sampled_inv, 0, 1);
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      std::cout << f_basis_sampled_inv[i][j] << " ";
    }
  }
  EXPECT_FALSE(0 == 1);
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
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
