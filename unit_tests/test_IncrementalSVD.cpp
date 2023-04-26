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
// Framework to run unit tests on the CAROM::IncrementalSVD class.

#include <iostream>

#ifdef CAROM_HAS_GTEST
#include<gtest/gtest.h>
#include <mpi.h>
#include "linalg/svd/IncrementalSVD.h"

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
    (
        CAROM::Options options,
        const std::string& basis_file_name)
        : CAROM::IncrementalSVD(
              options,
              basis_file_name)
    {
        int dim = options.dim;

        /* Construct a fake d_U, d_S, d_basis */
        d_basis = new CAROM::Matrix(dim, dim, false);
        d_S = new CAROM::Vector(dim, false);

        /* Use the identity matrix as a fake basis and fake singular values */
        for (int i = 0; i < dim; i++)
        {
            for (int j = 0; j < i; j++)
            {
                d_basis->item(i, j) = 0;
                d_basis->item(j, i) = 0;
            }
            d_basis->item(i, i) = d_S->item(i) = 1;
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
    CAROM::Options incremental_svd_options = CAROM::Options(3,
            4).setMaxBasisDimension(3)
            .setIncrementalSVD(1e-1, -1.0, -1.0, -1.0);

    FakeIncrementalSVD svd(
        incremental_svd_options,
        "irrelevant.txt");

    const CAROM::Matrix *B = svd.getSpatialBasis();
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
    CAROM::Options incremental_svd_options = CAROM::Options(3,
            4).setMaxBasisDimension(3)
            .setIncrementalSVD(1e-1, -1.0, -1.0, -1.0);

    FakeIncrementalSVD svd(
        incremental_svd_options,
        "irrelevant.txt");

    const CAROM::Vector *S = svd.getSingularValues();
    for (int i = 0; i < svd.getDim(); i++)
    {
        EXPECT_DOUBLE_EQ(S->item(i), 1);
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
