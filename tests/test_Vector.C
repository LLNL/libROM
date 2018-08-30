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
