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
#include "../SVD.h"

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
  SUCCEED();
}

/**
 *  Fake SVD to test parts of the SVD abstract base class. (For those
 *  unfamiliar, "fake" or "mock" is a term of art in unit testing
 *  referring to an implementation that simulates the behavior of
 *  real objects.)
 */
class FakeSVD : public CAROM::SVD
{
public:

  FakeSVD(int dim,
	  int samples_per_time_interval)
    : SVD(dim, samples_per_time_interval, false)
  {
  }

  ~FakeSVD()
  {
  }

  /**
   * Stub implementations of methods not really testable from the
   * abstract base class, because there is no meaningful data we could
   * return without implementing an actual singular value decomposition.
   *
   */
  virtual const CAROM::Matrix* getBasis() { return NULL; }
  virtual const CAROM::Matrix* getSingularValues() { return NULL; }

  /**
   *  The only testable methods from the SVD abstract base class are
   *  those with concrete implementations, namely:
   *
   *  int getDim() const;
   *
   *  int getNumBasisTimeIntervals() const;
   *
   *  double getBasisIntervalStartTime(int) const;
   *
   *  bool isNewTimeInterval() const;
   */

  /**
   * This method needs to do a few things:
   *
   * 1) If there are no sample time intervals stored, create one and
   *    record its start time. Set the number of samples in the current
   *    time interval to zero.
   *
   * 2) If adding another sample to the current sample time interval
   *    would exceed the number of samples per time interval set in
   *    the constructor, create a new time interval and record its
   *    start time. Set the number of samples in the current time interval
   *    to zero.
   *
   * 3) Increment the number of samples in the current sample time
   *    interval.
   *
   * Implementing this behavior suffices for testing the abstract base class.
   *
   */
  bool takeSample(__attribute__((unused)) double* u_in, double time)
  {
    /**
       If a new time interval is needed, add one and reset the number
       of samples counter to zero.
    */
    if (isNewTimeInterval())
    {
      d_time_interval_start_times.push_back(time);
      d_num_samples = 0;
    }

    /* Increment the number of samples in the current time interval */
    d_num_samples++;

    /**
	This method should almost always succeed because it does not
	do anything likely to fail, so it should return true.
    */
    return true;
  }

};

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
