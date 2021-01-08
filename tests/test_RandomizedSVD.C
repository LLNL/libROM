

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
// Framework to run unit tests on the CAROM::RandomizedSVD class.

#include <iostream>

#ifdef CAROM_HAS_GTEST
#include<gtest/gtest.h>
#include <mpi.h>
#include "../BasisGenerator.h"
#define _USE_MATH_DEFINES
#include <cmath>

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
  SUCCEED();
}

TEST(RandomizedSVDTest, Test_RandomizedSVD)
{
  // Get the rank of this process, and the number of processors.
  int mpi_init, d_rank, d_num_procs;
  MPI_Initialized(&mpi_init);
  if (mpi_init == 0) {
    MPI_Init(nullptr, nullptr);
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

  int num_total_rows = 5;
  int d_num_rows = num_total_rows / d_num_procs;
  if (num_total_rows % d_num_procs > d_rank) {
    d_num_rows++;
  }
  int *row_offset = new int[d_num_procs + 1];
  row_offset[d_num_procs] = num_total_rows;
  row_offset[d_rank] = d_num_rows;

  MPI_Allgather(MPI_IN_PLACE,
			     1,
			     MPI_INT,
			     row_offset,
			     1,
			     MPI_INT,
			     MPI_COMM_WORLD);

  for (int i = d_num_procs - 1; i >= 0; i--) {
    row_offset[i] = row_offset[i + 1] - row_offset[i];
  }

  double* sample1 = new double[5] {0.5377, 1.8339, -2.2588, 0.8622, 0.3188};
  double* sample2 = new double[5] {-1.3077, -0.4336, 0.3426, 3.5784, 2.7694};
  double* sample3 = new double[5] {-1.3499, 3.0349, 0.7254, -0.0631, 0.7147};

  double* basis_true_ans = new double[15] {
    -3.08158946098238906153e-01,     -9.49897947980622436859e-02,     4.50691774108525455844e-01,
    1.43697905723454588678e-01,      9.53289043424091042667e-01,      -8.77767692937215793236e-02,
    2.23655845793718777159e-02,      -2.10628953513210981363e-01,     -8.42235962392685610922e-01,
    7.29903965154318656872e-01,      -1.90917141788944810799e-01,     2.77280930877637610266e-01,
    5.92561353877168350834e-01,      -3.74570084880572404251e-02,     -5.40928141934187284301e-02};

   double* basis_right_true_ans = new double[9] {
     1.78651649346571128607e-01,      5.44387957786310439090e-01,      8.19588518467041615700e-01,
     9.49719639253861713790e-01,      -3.13100149275942318816e-01,     9.50441422536085767092e-04,
     2.57130696341889120049e-01,      7.78209514167381932737e-01,      -5.72951792961766015466e-01};

   double* sv_true_ans = new double[9] {
     4.84486375065219565528e+00,      0.00000000000000000000e+00,      0.00000000000000000000e+00,
     0.00000000000000000000e+00,      3.66719976398777092186e+00,      0.00000000000000000000e+00,
     0.00000000000000000000e+00,      0.00000000000000000000e+00,      2.69114625366671766926e+00};

  CAROM::Options randomized_svd_options = CAROM::Options(d_num_rows, 3, 1);
  randomized_svd_options.setMaxBasisDimension(num_total_rows);
  randomized_svd_options.setDebugMode(true);
  randomized_svd_options.setRandomizedSVD(true);
  CAROM::BasisGenerator sampler(randomized_svd_options, false);
  sampler.takeSample(&sample1[row_offset[d_rank]], 0, 0);
  sampler.takeSample(&sample2[row_offset[d_rank]], 0, 0);
  sampler.takeSample(&sample3[row_offset[d_rank]], 0, 0);

  const CAROM::Matrix* d_basis = sampler.getSpatialBasis();
  const CAROM::Matrix* d_basis_right = sampler.getTemporalBasis();
  const CAROM::Matrix* sv = sampler.getSingularValues();

  EXPECT_EQ(d_basis->numRows(), d_num_rows);
  EXPECT_EQ(d_basis->numColumns(), 3);
  EXPECT_EQ(d_basis_right->numRows(), 3);
  EXPECT_EQ(d_basis_right->numColumns(), 3);
  EXPECT_EQ(sv->numRows(), 3);
  EXPECT_EQ(sv->numColumns(), 3);

  double* d_basis_vals = d_basis->getData();
  double* d_basis_right_vals = d_basis_right->getData();
  double* sv_vals = sv->getData();

  for (int i = 0; i < d_num_rows * 3; i++) {

    EXPECT_NEAR(abs(d_basis_vals[i]), abs(basis_true_ans[row_offset[d_rank] * 3 + i]), 1e-7);
  }

  for (int i = 0; i < 9; i++) {
    EXPECT_NEAR(abs(d_basis_right_vals[i]), abs(basis_right_true_ans[i]), 1e-7);
  }

  for (int i = 0; i < 9; i++) {
    EXPECT_NEAR(sv_vals[i], sv_true_ans[i], 1e-7);
  }
}

TEST(RandomizedSVDTest, Test_RandomizedSVDTransposed)
{

  // Get the rank of this process, and the number of processors.
  int mpi_init, d_rank, d_num_procs;
  MPI_Initialized(&mpi_init);
  if (mpi_init == 0) {
    MPI_Init(nullptr, nullptr);
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

  int num_total_rows = 3;
  int d_num_rows = num_total_rows / d_num_procs;
  if (num_total_rows % d_num_procs > d_rank) {
    d_num_rows++;
  }
  int *row_offset = new int[d_num_procs + 1];
  row_offset[d_num_procs] = num_total_rows;
  row_offset[d_rank] = d_num_rows;

  MPI_Allgather(MPI_IN_PLACE,
           1,
           MPI_INT,
           row_offset,
           1,
           MPI_INT,
           MPI_COMM_WORLD);

  for (int i = d_num_procs - 1; i >= 0; i--) {
    row_offset[i] = row_offset[i + 1] - row_offset[i];
  }

  double* sample1 = new double[5] {0.5377, -1.3077, -1.3499};
  double* sample2 = new double[5] {1.8339, -0.4336, 3.0349};
  double* sample3 = new double[5] {-2.2588, 0.3426, 0.7254};
  double* sample4 = new double[5] {0.8622, 3.5784, -0.0631};
  double* sample5 = new double[5] {0.3188, 2.7694, 0.7147};

  double* basis_right_true_ans = new double[15] {
    -3.08158946098238906153e-01,     -9.49897947980622436859e-02,     4.50691774108525455844e-01,
    1.43697905723454588678e-01,      9.53289043424091042667e-01,      -8.77767692937215793236e-02,
    2.23655845793718777159e-02,      -2.10628953513210981363e-01,     -8.42235962392685610922e-01,
    7.29903965154318656872e-01,      -1.90917141788944810799e-01,     2.77280930877637610266e-01,
    5.92561353877168350834e-01,      -3.74570084880572404251e-02,     -5.40928141934187284301e-02};

   double* basis_true_ans = new double[9] {
     1.78651649346571128607e-01,      5.44387957786310439090e-01,      8.19588518467041615700e-01,
     9.49719639253861713790e-01,      -3.13100149275942318816e-01,     9.50441422536085767092e-04,
     2.57130696341889120049e-01,      7.78209514167381932737e-01,      -5.72951792961766015466e-01};

   double* sv_true_ans = new double[9] {
     4.84486375065219565528e+00,      0.00000000000000000000e+00,      0.00000000000000000000e+00,
     0.00000000000000000000e+00,      3.66719976398777092186e+00,      0.00000000000000000000e+00,
     0.00000000000000000000e+00,      0.00000000000000000000e+00,      2.69114625366671766926e+00};

  CAROM::Options randomized_svd_options = CAROM::Options(d_num_rows, 5, 1);
  randomized_svd_options.setMaxBasisDimension(num_total_rows);
  randomized_svd_options.setDebugMode(true);
  randomized_svd_options.setRandomizedSVD(true);
  CAROM::BasisGenerator sampler(randomized_svd_options, false);
  sampler.takeSample(&sample1[row_offset[d_rank]], 0, 0);
  sampler.takeSample(&sample2[row_offset[d_rank]], 0, 0);
  sampler.takeSample(&sample3[row_offset[d_rank]], 0, 0);
  sampler.takeSample(&sample4[row_offset[d_rank]], 0, 0);
  sampler.takeSample(&sample5[row_offset[d_rank]], 0, 0);

  const CAROM::Matrix* d_basis = sampler.getSpatialBasis();
  const CAROM::Matrix* d_basis_right = sampler.getTemporalBasis();
  const CAROM::Matrix* sv = sampler.getSingularValues();

  num_total_rows = 5;
  d_num_rows = num_total_rows / d_num_procs;
  if (num_total_rows % d_num_procs > d_rank) {
    d_num_rows++;
  }

  row_offset[d_num_procs] = num_total_rows;
  row_offset[d_rank] = d_num_rows;

  EXPECT_EQ(d_basis_right->numRows(), d_num_rows);
  EXPECT_EQ(d_basis_right->numColumns(), 3);
  EXPECT_EQ(d_basis->numRows(), 3);
  EXPECT_EQ(d_basis->numColumns(), 3);
  EXPECT_EQ(sv->numRows(), 3);
  EXPECT_EQ(sv->numColumns(), 3);

  double* d_basis_vals = d_basis->getData();
  double* d_basis_right_vals = d_basis_right->getData();
  double* sv_vals = sv->getData();

  MPI_Allgather(MPI_IN_PLACE,
			     1,
			     MPI_INT,
			     row_offset,
			     1,
			     MPI_INT,
			     MPI_COMM_WORLD);

  for (int i = d_num_procs - 1; i >= 0; i--) {
    row_offset[i] = row_offset[i + 1] - row_offset[i];
  }

  for (int i = 0; i < d_num_rows * 3; i++) {
    EXPECT_NEAR(abs(d_basis_right_vals[i]), abs(basis_right_true_ans[row_offset[d_rank] * 3 + i]), 1e-7);
  }

  for (int i = 0; i < 9; i++) {
    EXPECT_NEAR(abs(d_basis_vals[i]), abs(basis_true_ans[i]), 1e-7);
  }

  for (int i = 0; i < 9; i++) {
    EXPECT_NEAR(sv_vals[i], sv_true_ans[i], 1e-7);
  }
}

TEST(RandomizedSVDTest, Test_RandomizedSVDSmallerSubspace)
{

  // Get the rank of this process, and the number of processors.
  int mpi_init, d_rank, d_num_procs;
  MPI_Initialized(&mpi_init);
  if (mpi_init == 0) {
    MPI_Init(nullptr, nullptr);
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

  int num_total_rows = 5;
  int d_num_rows = num_total_rows / d_num_procs;
  if (num_total_rows % d_num_procs > d_rank) {
    d_num_rows++;
  }
  int *row_offset = new int[d_num_procs + 1];
  row_offset[d_num_procs] = num_total_rows;
  row_offset[d_rank] = d_num_rows;

  MPI_Allgather(MPI_IN_PLACE,
			     1,
			     MPI_INT,
			     row_offset,
			     1,
			     MPI_INT,
			     MPI_COMM_WORLD);

  for (int i = d_num_procs - 1; i >= 0; i--) {
    row_offset[i] = row_offset[i + 1] - row_offset[i];
  }

  double* sample1 = new double[5] {0.5377, 1.8339, -2.2588, 0.8622, 0.3188};
  double* sample2 = new double[5] {-1.3077, -0.4336, 0.3426, 3.5784, 2.7694};
  double* sample3 = new double[5] {-1.3499, 3.0349, 0.7254, -0.0631, 0.7147};

double* basis_true_ans = new double[10] {
  -3.66686936855722067641e-01,     -3.66854732104725300701e-01,
  -1.03335953010555314791e-01,     8.34110828036275564479e-01,
  2.73086819017167137247e-01,     3.22486279456235425123e-01,
  6.69406991793326300311e-01,      -2.50470054687995880016e-01,
  5.76350795759026368614e-01,      5.42594079591339900626e-02};

  double* basis_right_true_ans = new double[4] {
    -5.18839658032825834511e-02,     1.19929948552507714687e-01,
    9.85719269437193235106e-01,      -1.53109727976926013326e-01};

 double* sv_true_ans = new double[4] {
   4.67621414826046599700e+00,      0.00000000000000000000e+00,
   0.00000000000000000000e+00,      3.37971184883141617661e+00};

  CAROM::Options randomized_svd_options = CAROM::Options(d_num_rows, 3, 1);
  randomized_svd_options.setMaxBasisDimension(num_total_rows);
  randomized_svd_options.setDebugMode(true);
  randomized_svd_options.setRandomizedSVD(true, 2);
  CAROM::BasisGenerator sampler(randomized_svd_options, false);
  sampler.takeSample(&sample1[row_offset[d_rank]], 0, 0);
  sampler.takeSample(&sample2[row_offset[d_rank]], 0, 0);
  sampler.takeSample(&sample3[row_offset[d_rank]], 0, 0);

  const CAROM::Matrix* d_basis = sampler.getSpatialBasis();
  const CAROM::Matrix* d_basis_right = sampler.getTemporalBasis();
  const CAROM::Matrix* sv = sampler.getSingularValues();

  EXPECT_EQ(d_basis->numRows(), d_num_rows);
  EXPECT_EQ(d_basis->numColumns(), 2);
  EXPECT_EQ(d_basis_right->numRows(), 2);
  EXPECT_EQ(d_basis_right->numColumns(), 2);
  EXPECT_EQ(sv->numRows(), 2);
  EXPECT_EQ(sv->numColumns(), 2);

  double* d_basis_vals = d_basis->getData();
  double* d_basis_right_vals = d_basis_right->getData();
  double* sv_vals = sv->getData();

  for (int i = 0; i < d_num_rows * 2; i++) {
    EXPECT_NEAR(abs(d_basis_vals[i]), abs(basis_true_ans[row_offset[d_rank] * 2 + i]), 1e-7);
  }

  for (int i = 0; i < 4; i++) {
    EXPECT_NEAR(abs(d_basis_right_vals[i]), abs(basis_right_true_ans[i]), 1e-7);
  }

  for (int i = 0; i < 4; i++) {
    EXPECT_NEAR(sv_vals[i], sv_true_ans[i], 1e-7);
  }
}

TEST(RandomizedSVDTest, Test_RandomizedSVDTransposedSmallerSubspace)
{

  // Get the rank of this process, and the number of processors.
  int mpi_init, d_rank, d_num_procs;
  MPI_Initialized(&mpi_init);
  if (mpi_init == 0) {
    MPI_Init(nullptr, nullptr);
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

  int num_total_rows = 3;
  int d_num_rows = num_total_rows / d_num_procs;
  if (num_total_rows % d_num_procs > d_rank) {
    d_num_rows++;
  }
  int *row_offset = new int[d_num_procs + 1];
  row_offset[d_num_procs] = num_total_rows;
  row_offset[d_rank] = d_num_rows;

  MPI_Allgather(MPI_IN_PLACE,
           1,
           MPI_INT,
           row_offset,
           1,
           MPI_INT,
           MPI_COMM_WORLD);

  for (int i = d_num_procs - 1; i >= 0; i--) {
    row_offset[i] = row_offset[i + 1] - row_offset[i];
  }

  double* sample1 = new double[5] {0.5377, -1.3077, -1.3499};
  double* sample2 = new double[5] {1.8339, -0.4336, 3.0349};
  double* sample3 = new double[5] {-2.2588, 0.3426, 0.7254};
  double* sample4 = new double[5] {0.8622, 3.5784, -0.0631};
  double* sample5 = new double[5] {0.3188, 2.7694, 0.7147};

  double* basis_right_true_ans = new double[10] {
    -3.66686936855722067641e-01,     -3.66854732104725300701e-01,
    -1.03335953010555314791e-01,     8.34110828036275564479e-01,
    2.73086819017167137247e-01,     3.22486279456235425123e-01,
    6.69406991793326300311e-01,      -2.50470054687995880016e-01,
    5.76350795759026368614e-01,      5.42594079591339900626e-02};

   double* basis_true_ans = new double[4] {
      -5.18839658032825834511e-02,     1.19929948552507714687e-01,
      9.85719269437193235106e-01,      -1.53109727976926013326e-01};

   double* sv_true_ans = new double[4] {
     4.67621414826046599700e+00,      0.00000000000000000000e+00,
     0.00000000000000000000e+00,      3.37971184883141617661e+00};

  CAROM::Options randomized_svd_options = CAROM::Options(d_num_rows, 5, 1);
  randomized_svd_options.setMaxBasisDimension(num_total_rows);
  randomized_svd_options.setDebugMode(true);
  randomized_svd_options.setRandomizedSVD(true, 2);
  CAROM::BasisGenerator sampler(randomized_svd_options, false);
  sampler.takeSample(&sample1[row_offset[d_rank]], 0, 0);
  sampler.takeSample(&sample2[row_offset[d_rank]], 0, 0);
  sampler.takeSample(&sample3[row_offset[d_rank]], 0, 0);
  sampler.takeSample(&sample4[row_offset[d_rank]], 0, 0);
  sampler.takeSample(&sample5[row_offset[d_rank]], 0, 0);

  const CAROM::Matrix* d_basis = sampler.getSpatialBasis();
  const CAROM::Matrix* d_basis_right = sampler.getTemporalBasis();
  const CAROM::Matrix* sv = sampler.getSingularValues();

  num_total_rows = 5;
  d_num_rows = num_total_rows / d_num_procs;
  if (num_total_rows % d_num_procs > d_rank) {
    d_num_rows++;
  }

  row_offset[d_num_procs] = num_total_rows;
  row_offset[d_rank] = d_num_rows;

  EXPECT_EQ(d_basis_right->numRows(), d_num_rows);
  EXPECT_EQ(d_basis_right->numColumns(), 2);
  EXPECT_EQ(d_basis->numRows(), 2);
  EXPECT_EQ(d_basis->numColumns(), 2);
  EXPECT_EQ(sv->numRows(), 2);
  EXPECT_EQ(sv->numColumns(), 2);

  double* d_basis_vals = d_basis->getData();
  double* d_basis_right_vals = d_basis_right->getData();
  double* sv_vals = sv->getData();

  MPI_Allgather(MPI_IN_PLACE,
			     1,
			     MPI_INT,
			     row_offset,
			     1,
			     MPI_INT,
			     MPI_COMM_WORLD);

  for (int i = d_num_procs - 1; i >= 0; i--) {
    row_offset[i] = row_offset[i + 1] - row_offset[i];
  }

  for (int i = 0; i < d_num_rows * 2; i++) {
    EXPECT_NEAR(abs(d_basis_right_vals[i]), abs(basis_right_true_ans[row_offset[d_rank] * 2 + i]), 1e-7);
  }

  for (int i = 0; i < 4; i++) {
    EXPECT_NEAR(abs(d_basis_vals[i]), abs(basis_true_ans[i]), 1e-7);
  }

  for (int i = 0; i < 4; i++) {
    EXPECT_NEAR(sv_vals[i], sv_true_ans[i], 1e-7);
  }
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
