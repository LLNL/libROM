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
#include "../Matrix.h"

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
  SUCCEED();
}

TEST(DEIMSerialTest, Test_DEIM)
{

  // Orthonormal input matrix to DEIM
  double* orthonormal_mat = new double[25] {
   -0.1735,   -0.2978,    0.4339,    0.8304,    0.0575,
   -0.5919,   -0.1608,   -0.7489,    0.2185,   -0.1230,
    0.7291,    0.1574,   -0.4449,    0.4549,   -0.1970,
   -0.2783,    0.7273,    0.1992,    0.1386,   -0.5785,
   -0.1029,    0.5760,   -0.1147,    0.1910,    0.7798};

   // Result of DEIM (f_basis_sampled_inv)
   double* DEIM_true_ans = new double[25] {
  -0.173528,  -0.591896,  0.729022,  -0.278257,  -0.102821,
  -0.29779,   -0.16083,   0.157456,   0.727281 ,  0.575906,
  0.433945,   -0.748987, -0.444923,   0.19924,   -0.11473,
  0.830464,   0.218532,   0.454906,   0.138637,   0.191005,
  0.0574838, -0.123002,  -0.196971,  -0.578576,   0.779759};

  CAROM::Matrix* u = new CAROM::Matrix(orthonormal_mat, 5, 5, false);
  double* DEIM_res = NULL;
  int* f_sampled_row = new int[5] {0};
  int* f_sampled_rows_per_proc = new int[5] {0};
  CAROM::Matrix f_basis_sampled_inv = CAROM::Matrix(5, 5, false);
  CAROM::DEIM(u, 5, f_sampled_row, f_sampled_rows_per_proc, f_basis_sampled_inv, 0, 1);

  // Compare the norm between the DEIM result and the true DEIM answer
  double l2_norm_diff = 0.0;
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      l2_norm_diff += pow(abs(DEIM_true_ans[i * 5 + j] - f_basis_sampled_inv(i, j)), 2);
    }
  }
  l2_norm_diff = sqrt(l2_norm_diff);

  // Allow for some error due to float rounding
  EXPECT_TRUE(l2_norm_diff < 1e-5);
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
