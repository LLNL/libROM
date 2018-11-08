/******************************************************************************
 *
 * Copyright (c) 2013-2018, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 * Written by William Arrighi wjarrighi@llnl.gov
 *            Geoff Oxberry oxberry1@llnl.gov
 * CODE-686965
 * All rights reserved.
 *
 * This file is part of libROM.
 * For details, see https://computation.llnl.gov/librom
 * Please also read README_BSD_NOTICE.
 *
 * Redistribution and use in source and binary forms, with or without
 * modifications, are permitted provided that the following conditions are met:
 *
 *    o Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the disclaimer below.
 *    o Redistribution in binary form must reproduce the above copyright
 *      notice, this list of conditions and the disclaimer (as noted below) in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *    o Neither the name of the LLNS/LLNL nor the names of its contributors may
 *      be used to endorse or promote products derived from this software
 *      without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
 * LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OR SUCH DAMAGE.
 *
 *****************************************************************************/

#include <iostream>

#ifdef CAROM_HAS_GTEST
#include<gtest/gtest.h>
#include <mpi.h>
#include "../IncrementalSVD.h"

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
  SUCCEED();
}

/**
 *  Fake/mock version of an IncrementalSVD. The only public methods
 *  to test are:
 *
 * * getBasis()
 *
 * * getSingularValues()
 *
 */
class FakeIncrementalSVD : public CAROM::IncrementalSVD
{
public:

  FakeIncrementalSVD
  (int dim,
   double linearity_tol,
   bool skip_linearly_dependent,
   int max_basis_dimension,
   int samples_per_time_interval,
   const std::string& basis_file_name)
    : CAROM::IncrementalSVD(dim,
			    linearity_tol,
			    skip_linearly_dependent,
			    max_basis_dimension,
			    samples_per_time_interval,
			    basis_file_name,
			    false,
			    false,
			    false)
  {
    /* Construct a fake d_U, d_S, d_basis */
    d_basis = new CAROM::Matrix(dim, dim, false);
    d_S = new CAROM::Matrix(dim, dim, false);

    /* Use the identity matrix as a fake basis and fake singular values */
    for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < i; j++)
      {
	d_basis->item(i, j) = d_S->item(i, j) = 0;
	d_basis->item(j, i) = d_S->item(j, i) = 0;
      }
      d_basis->item(i, i) = d_S->item(i, i) = 1;
    }
  }

  ~FakeIncrementalSVD()
  {
  }

  void buildInitialSVD
  (__attribute__((unused)) double* u,
   __attribute__((unused)) double time)
  {
    /* Do nothing */
  }

  void computeBasis()
  {
    /* Do nothing */
  }

  void addLinearlyDependentSample
  (__attribute__((unused)) const CAROM::Matrix *A,
   __attribute__((unused)) const CAROM::Matrix *W,
   __attribute__((unused)) const CAROM::Matrix *sigma)
  {
    /* Do nothing */
  }

  void addNewSample
  (__attribute__((unused)) const CAROM::Vector *j,
   __attribute__((unused)) const CAROM::Matrix *A,
   __attribute__((unused)) const CAROM::Matrix *W,
   __attribute__((unused)) CAROM::Matrix *sigma)
  {
    /* Do nothing */
  }

};

TEST(IncrementalSVDSerialTest, Test_getBasis)
{
  FakeIncrementalSVD svd(3,
			 1e-1,
			 false,
			 3,
			 4,
			 "irrelevant.txt");

  const CAROM::Matrix *B = svd.getBasis();
  for (int i = 0; i < svd.getDim(); i++)
  {
    for (int j = 0; j < i; j++)
    {
      EXPECT_DOUBLE_EQ(B->item(i, j), 0);
      EXPECT_DOUBLE_EQ(B->item(j, i), 0);
    }
    EXPECT_DOUBLE_EQ(B->item(i, i), 1);
  }
}

TEST(IncrementalSVDSerialTest, Test_getSingularValues)
{
  FakeIncrementalSVD svd(3,
			 1e-1,
			 false,
			 3,
			 4,
			 "irrelevant.txt");

  const CAROM::Matrix *S = svd.getSingularValues();
  for (int i = 0; i < svd.getDim(); i++)
  {
    for (int j = 0; j < i; j++)
    {
      EXPECT_DOUBLE_EQ(S->item(i, j), 0);
      EXPECT_DOUBLE_EQ(S->item(j, i), 0);
    }
    EXPECT_DOUBLE_EQ(S->item(i, i), 1);
  }
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
