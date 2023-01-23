/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: This source file is a test runner that uses the Google Test
// Framework to run unit tests on the CAROM::SVD class.

#include <iostream>

#ifdef CAROM_HAS_GTEST
#include<gtest/gtest.h>
#include <mpi.h>
#include "linalg/Options.h"
#include "linalg/svd/SVD.h"

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

    FakeSVD(CAROM::Options options)
        : SVD(options)
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
    virtual const CAROM::Matrix* getSpatialBasis() {
        return NULL;
    }
    virtual const CAROM::Matrix* getTemporalBasis() {
        return NULL;
    }
    virtual const CAROM::Vector* getSingularValues() {
        return NULL;
    }
    virtual const CAROM::Matrix* getSnapshotMatrix() {
        return NULL;
    }

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
    bool takeSample
    (__attribute__((unused)) double* u_in,
     double time,
     __attribute__((unused)) bool add_without_increase)
    {
        /**
           If a new time interval is needed, add one and reset the number
           of samples counter to zero.
        */
        if (isNewTimeInterval())
        {
            int num_time_intervals =
                static_cast<int>(d_time_interval_start_times.size());
            increaseTimeInterval();
            d_time_interval_start_times[num_time_intervals] = time;
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

TEST(SVDSerialTest, Test_getDim)
{
    FakeSVD svd(CAROM::Options(5, 2));
    EXPECT_EQ(svd.getDim(), 5);
}

TEST(SVDSerialTest, Test_isNewTimeInterval)
{
    FakeSVD svd(CAROM::Options(5, 2));

    /* 0 samples, so taking a sample will create a new time interval */
    EXPECT_TRUE(svd.isNewTimeInterval());

    /* 1 sample; limit is 2, taking a sample won't create a new time interval */
    svd.takeSample(NULL, 0, true);
    EXPECT_FALSE(svd.isNewTimeInterval());

    /* 2 samples; limit is 2, taking a sample will create a new time interval */
    svd.takeSample(NULL, 0.5, true);
    EXPECT_TRUE(svd.isNewTimeInterval());

    /* 1 sample; limit is 2, taking a sample won't create a new time interval */
    svd.takeSample(NULL, 1, true);
    EXPECT_FALSE(svd.isNewTimeInterval());

    /* 2 samples; limit is 2, taking a sample will create a new time interval */
    svd.takeSample(NULL, 1.5, true);
    EXPECT_TRUE(svd.isNewTimeInterval());

    /* 1 sample; limit is 2, taking a sample won't create a new time interval */
    svd.takeSample(NULL, 2, true);
    EXPECT_FALSE(svd.isNewTimeInterval());
}

TEST(SVDSerialTest, Test_getNumBasisTimeIntervals)
{
    FakeSVD svd(CAROM::Options(5, 2));

    /* Number of time intervals starts at zero. */
    EXPECT_EQ(svd.getNumBasisTimeIntervals(), 0);

    /* Creates new time interval; number of intervals = 1 */
    svd.takeSample(NULL, 0, true);
    EXPECT_EQ(svd.getNumBasisTimeIntervals(), 1);
    svd.takeSample(NULL, 0.5, true);
    EXPECT_EQ(svd.getNumBasisTimeIntervals(), 1);

    /* Creates new time interval; number of intervals = 2 */
    svd.takeSample(NULL, 1, true);
    EXPECT_EQ(svd.getNumBasisTimeIntervals(), 2);
    svd.takeSample(NULL, 1.5, true);
    EXPECT_EQ(svd.getNumBasisTimeIntervals(), 2);

    /* Creates new time interval; number of intervals = 3 */
    svd.takeSample(NULL, 2, true);
    EXPECT_EQ(svd.getNumBasisTimeIntervals(), 3);
}

TEST(SVDSerialTest, Test_getBasisIntervalStartTime)
{
    FakeSVD svd(CAROM::Options(5, 2));

    /* 1st time interval starts at time 0 */
    svd.takeSample(NULL, 0, true);
    EXPECT_DOUBLE_EQ(svd.getBasisIntervalStartTime(0), 0);

    svd.takeSample(NULL, 0.5, true);
    EXPECT_DOUBLE_EQ(svd.getBasisIntervalStartTime(0), 0);

    /* 2nd time interval starts at time 1 */
    svd.takeSample(NULL, 1, true);
    EXPECT_DOUBLE_EQ(svd.getBasisIntervalStartTime(0), 0);
    EXPECT_DOUBLE_EQ(svd.getBasisIntervalStartTime(1), 1);

    svd.takeSample(NULL, 1.5, true);
    EXPECT_DOUBLE_EQ(svd.getBasisIntervalStartTime(0), 0);
    EXPECT_DOUBLE_EQ(svd.getBasisIntervalStartTime(1), 1);

    /* 3rd time interval starts at time 2 */
    svd.takeSample(NULL, 2, true);
    EXPECT_DOUBLE_EQ(svd.getBasisIntervalStartTime(0), 0);
    EXPECT_DOUBLE_EQ(svd.getBasisIntervalStartTime(1), 1);
    EXPECT_DOUBLE_EQ(svd.getBasisIntervalStartTime(2), 2);
}


TEST(SVDSerialTest, Test_increaseTimeInterval)
{
    FakeSVD svd(CAROM::Options(5, 2, 2));

    ASSERT_NO_THROW(svd.takeSample(NULL, 0, true));
    ASSERT_NO_THROW(svd.takeSample(NULL, 0.5, true));
    ASSERT_NO_THROW(svd.takeSample(NULL, 1, true));
    ASSERT_NO_THROW(svd.takeSample(NULL, 1.5, true));

    /* The maximum number of time intervals is surpassed */
    EXPECT_DEATH(svd.takeSample(NULL, 2, true), ".*");
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
