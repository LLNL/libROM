/******************************************************************************
 *
 * Copyright (c) 2013-2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 * Written by William Arrighi wjarrighi@llnl.gov
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
#include "../Vector.h"
#define _USE_MATH_DEFINES
#include <cmath>

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
  SUCCEED();
}

/**
 *  Test methods that do not require assigning data to the Vector. These
 *  methods are:
 *
 *  TODO(oxberry1@llnl.gov): Do more exhaustive testing to test all branches
 *  of each method. For now, simple tests are better than nothing.
 *
 *  * Vector(int, bool)
 *
 *  * bool distributed() const
 *
 *  * int dim() const
 *
 *  * void setSize(int)
 *
 */

TEST(VectorSerialTest, Test_distributed)
{
  CAROM::Vector v(2, false);
  EXPECT_FALSE(v.distributed());

  CAROM::Vector w(2, true);
  EXPECT_TRUE(w.distributed());
}

TEST(VectorSerialTest, Test_dim)
{
  CAROM::Vector v(2, false);
  EXPECT_EQ(v.dim(), 2);
}

TEST(VectorSerialTest, Test_setSize)
{
  CAROM::Vector v(2, false);
  EXPECT_EQ(v.dim(), 2);
  v.setSize(3);
  EXPECT_EQ(v.dim(), 3);
}

/** Test methods that require assigning data to the Vector
 *
 *  TODO(oxberry1@llnl.gov): Do more exhaustive testing to test all branches
 *  of each method. For now, simple tests are better than nothing.
 *
 *  * Vector(double*, int, bool, bool)
 *
 *  * const double& operator() (int, int) const
 *
 *  * double& operator() (int, int)
 *
 *  * const double& item(int, int) const
 *
 *  * double& item(int, int)
 *
 *  * Vector(const Vector&)
 *
 *  * Vector& operator=(const Vector&) for copy-assignment
 *
 *  * Vector& operator=(const Vector&) for assignment
 */

TEST(VectorSerialTest, Test_const_call_operator)
{
  double v_data[2] = {1, 2};
  const CAROM::Vector v(v_data, 2, false, true);

  EXPECT_DOUBLE_EQ(v(0), 1);
  EXPECT_DOUBLE_EQ(v(1), 2);
}

TEST(VectorSerialTest, Test_nonconst_call_operator)
{
  CAROM::Vector v(2, false);
  v(0) = 1;
  v(1) = 2;

  EXPECT_DOUBLE_EQ(v(0), 1);
  EXPECT_DOUBLE_EQ(v(1), 2);
}

TEST(VectorSerialTest, Test_const_item)
{
  double v_data[2] = {1, 2};
  const CAROM::Vector v(v_data, 2, false, true);

  EXPECT_DOUBLE_EQ(v.item(0), 1);
  EXPECT_DOUBLE_EQ(v.item(1), 2);
}

TEST(VectorSerialTest, Test_nonconst_item)
{
  CAROM::Vector v(2, false);
  v.item(0) = 1;
  v.item(1) = 2;

  EXPECT_DOUBLE_EQ(v.item(0), 1);
  EXPECT_DOUBLE_EQ(v.item(1), 2);
}

TEST(VectorSerialTest, Test_copy_constructor)
{
  double v_data[2] = {1, 2};
  const CAROM::Vector v(v_data, 2, false, true);
  CAROM::Vector w(v);

  EXPECT_FALSE(w.distributed());
  EXPECT_EQ(w.dim(), 2);
  EXPECT_DOUBLE_EQ(w(0), 1);
  EXPECT_DOUBLE_EQ(w(1), 2);
}

TEST(VectorSerialTest, Test_copy_assignment_operator)
{
  double v_data[2] = {1, 2};
  const CAROM::Vector v(v_data, 2, false, true);
  CAROM::Vector w = v;

  EXPECT_FALSE(w.distributed());
  EXPECT_EQ(w.dim(), 2);
  EXPECT_DOUBLE_EQ(w(0), 1);
  EXPECT_DOUBLE_EQ(w(1), 2);
}

TEST(VectorSerialTest, Test_assignment_operator)
{
  double v_data[2] = {1, 2};
  const CAROM::Vector v(v_data, 2, false, true);
  CAROM::Vector w(2, false);
  w = v;

  EXPECT_FALSE(w.distributed());
  EXPECT_EQ(w.dim(), 2);
  EXPECT_DOUBLE_EQ(w(0), 1);
  EXPECT_DOUBLE_EQ(w(1), 2);
}

/** Test methods that operate on Vector objects
 *
 *  * double norm() const
 *
 *  * double normalize()
 *
 *  * Vector* plus(const Vector&) const
 *
 *  * Vector* plus(const Vector*) const
 *
 *  * void plus(const Vector&, Vector*&) const
 *
 *  * void plus(const Vector&, Vector&) const
 *
 *  * Vector* plusAx(double, const Vector&)
 *
 *  * Vector* plusAx(double, const Vector*)
 *
 *  * void plusAx(double, const Vector&, Vector*&) const
 *
 *  * void plusAx(double, const Vector&, Vector&) const
 *
 *  * void plusEqAx(double, const Vector&)
 *
 *  * void plusEqAx(double, const Vector*)
 *
 *  * Vector* minus(const Vector&) const
 *
 *  * Vector* minus(const Vector*) const
 *
 *  * void minus(const Vector&, Vector*&) const
 *
 *  * void minus(const Vector&, Vector&) const
 *
 *  * Vector* mult(double) const
 *
 *  * void mult(double, Vector*&) const
 *
 *  * void mult(double, Vector&) const
 */

TEST(VectorSerialTest, TestNorm)
{
  CAROM::Vector v(2, false);

  /* 2-norm of (1, 1) = sqrt(2) */
  v(0) = 1; v(1) = 1;
  EXPECT_DOUBLE_EQ(v.norm(), M_SQRT2);

  /* 2-norm of (-1, 1) = sqrt(2) */
  v(0) = -1; v(1) = 1;
  EXPECT_DOUBLE_EQ(v.norm(), M_SQRT2);

  /* 2-norm of (3, 4) = 5 */
  v(0) = 3; v(1) = 4;
  EXPECT_DOUBLE_EQ(v.norm(), 5);

  /* 2-norm of (5, 12) = 13 */
  v(0) = 5; v(1) = 12;
  EXPECT_DOUBLE_EQ(v.norm(), 13);

  /** TODO(oxberry1@llnl.gov): Add more thorough tests **/
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
