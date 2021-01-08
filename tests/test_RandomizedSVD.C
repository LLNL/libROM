

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
    4.61558975111678149172e-01,      4.37365354601744860119e-01,      1.66627712118903870131e-01,
    -3.25109650190036403306e-01,     7.57728478636118030209e-01,      4.51858262649965980451e-01,
    -7.38549782321004655294e-03,     -4.81769074950681830405e-01,     8.71483680739033150786e-01,
    -6.23073359704019225092e-01,     -3.74247217193117640011e-02,     -6.99081595951640683007e-02,
    -5.41287419672967340389e-01,     -3.25122728101360919384e-02,     -6.07318652483036999778e-02};

   double* basis_right_true_ans = new double[9] {
     -2.19374126997083473967e-01,     8.20726595609202469461e-01,      -5.27525210453487769513e-01,
     -8.83662121247607013075e-01,     -3.96326568787719990539e-01,     -2.49131504120077035269e-01,
     -4.13541107843521604792e-01,     4.11501040257107097986e-01,      8.12188799473910205684e-01};

   double* sv_true_ans = new double[9] {
     4.74592085430968513293e+00,      0.00000000000000000000e+00,      0.00000000000000000000e+00,
     0.00000000000000000000e+00,      3.25364999902110074714e+00,      0.00000000000000000000e+00,
     0.00000000000000000000e+00,      0.00000000000000000000e+00,      2.14185949946548248590e+00};

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
    4.61558975111678149172e-01,      4.37365354601744860119e-01,      1.66627712118903870131e-01,
    -3.25109650190036403306e-01,     7.57728478636118030209e-01,      4.51858262649965980451e-01,
    -7.38549782321004655294e-03,     -4.81769074950681830405e-01,     8.71483680739033150786e-01,
    -6.23073359704019225092e-01,     -3.74247217193117640011e-02,     -6.99081595951640683007e-02,
    -5.41287419672967340389e-01,     -3.25122728101360919384e-02,     -6.07318652483036999778e-02};

   double* basis_true_ans = new double[9] {
     -2.19374126997083473967e-01,     8.20726595609202469461e-01,      -5.27525210453487769513e-01,
     -8.83662121247607013075e-01,     -3.96326568787719990539e-01,     -2.49131504120077035269e-01,
     -4.13541107843521604792e-01,     4.11501040257107097986e-01,      8.12188799473910205684e-01};

   double* sv_true_ans = new double[9] {
     4.74592085430968513293e+00,      0.00000000000000000000e+00,      0.00000000000000000000e+00,
     0.00000000000000000000e+00,      3.25364999902110074714e+00,      0.00000000000000000000e+00,
     0.00000000000000000000e+00,      0.00000000000000000000e+00,      2.14185949946548248590e+00};

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
    4.24468361793615245325e-01,      4.77458702069573692750e-01,
    -3.48848060484092337141e-01,     8.70309584749493048150e-01,
    1.68077499423117754374e-01,      2.43015657433558601719e-02,
    -6.17869712566928130926e-01,     -8.93349882780821724637e-02,
    -5.36766814373677014771e-01,     -7.76086869041048010853e-02};

   double* basis_right_true_ans = new double[4] {
     -3.18552083217684134375e-01,     5.63151606900417212032e-01,
     -8.61623137159813645702e-01,     -5.07334680179376329434e-01};

   double* sv_true_ans = new double[4] {
      4.69316598773643711695e+00,      0.00000000000000000000e+00,
      0.00000000000000000000e+00,      3.01185616415252432887e+00};

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
    4.24468361793615245325e-01,      4.77458702069573692750e-01,
    -3.48848060484092337141e-01,     8.70309584749493048150e-01,
    1.68077499423117754374e-01,      2.43015657433558601719e-02,
    -6.17869712566928130926e-01,     -8.93349882780821724637e-02,
    -5.36766814373677014771e-01,     -7.76086869041048010853e-02};

   double* basis_true_ans = new double[4] {
     -3.18552083217684134375e-01,     5.63151606900417212032e-01,
     -8.61623137159813645702e-01,     -5.07334680179376329434e-01};

   double* sv_true_ans = new double[4] {
      4.69316598773643711695e+00,      0.00000000000000000000e+00,
      0.00000000000000000000e+00,      3.01185616415252432887e+00};

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
