/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: This source file is a test runner that uses the Google Test
// Framework to run unit tests on the CAROM::Matrix class.

#include <iostream>
#include <cmath>

#ifdef CAROM_HAS_GTEST
#include<gtest/gtest.h>
#include "algo/manifold_interp/PCHIPInterpolator.h"
#include "linalg/Vector.h"
#include <cfloat>

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
    SUCCEED();
}

TEST(InterpolationTest,test_accuracy)
{
    int n_out = 25;
    int n_snap = 11;
    bool SUCCESS = true;
    //input times
    std::vector<double> t{ 0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15};
    //Test function from original PCHIP paper
    std::vector<double> y{10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85};
    //sin(x)
    std::vector<double> y2{0, 0.909297426825682, 0.141120008059867, -0.958924274663138, -0.279415498198926,
                           0.989358246623382, 0.412118485241757, -0.999990206550703, -0.536572918000435,
                           0.990607355694870, 0.650287840157117};

    //query points
    std::vector<double> tq{0, 0.6000, 1.2000, 1.8000, 2.4000, 3.0000, 3.6000, 4.2000, 4.8000, 5.4000, 6.0000, 6.6000,
                           7.2000, 7.8000, 8.4000, 9.0000, 9.6000, 10.2000, 10.8000, 11.4000, 12.0000, 12.6000,
                           13.2000, 13.8000, 14.4000, 15.000};
    //pchip reference solution for original test function.
    std::vector<double> yq1{10.000000000000000, 10.000000000000000, 10.000000000000002, 10.000000000000004,
                            10.000000000000000, 10.000000000000000, 10.000000000000002, 10.000000000000000,
                            10.000000000000004, 10.000000000000000, 10.000000000000000, 10.000000000000002,
                            10.000000000000000, 10.000000000000004, 10.102641509433962, 10.500000000000000,
                            11.106230625292378, 12.213163262123807, 14.128630750038976, 27.078413223140512,
                            50.000000000000000, 53.832363636363638, 55.720727272727267, 58.433818181818161,
                            67.055999999999969, 85.000000000000000};
    //pchip reference solution for sin(x)
    std::vector<double> yq2{0, 0.569748887844739, 0.833039030477175, 0.906694689302138, 0.701592409322499, 0.141120008059867,
                            -0.288488198334352, -0.697095554949408, -0.939878053603591, -0.782971083698006, -0.279415498198926,
                            0.188293444380395, 0.669217685146469, 0.965688937709033, 0.846474732609843, 0.412118485241757,
                            -0.077580693288345, -0.623537711025344, -0.971758328554163, -0.890773577229575, -0.536572918000435,
                            -0.041614069121016, 0.560852411851254, 0.957953731078007, 0.938810668593539, 0.650287840157117};

    std::vector<std::shared_ptr<CAROM::Vector>> snapshots;
    std::vector<std::shared_ptr<CAROM::Vector>> out_snapshots;
    std::vector<CAROM::Vector> reference_snapshots;
    std::vector<CAROM::Vector> times;
    std::vector<CAROM::Vector> out_times;
    for(int i = 0; i < t.size(); ++i)
    {
        CAROM::Vector temp_t(1, false);
        temp_t(0) = t[i];
        times.push_back(temp_t);
        CAROM::Vector *temp_y = new CAROM::Vector(2,false);
        (*temp_y)(0) = y[i];
        (*temp_y)(1) = y2[i];
        snapshots.push_back(std::shared_ptr<CAROM::Vector>(temp_y));
    }

    for(int i = 0; i < tq.size(); ++i)
    {
        CAROM::Vector temp_t(1, false);
        temp_t(0) = tq[i];
        out_times.push_back(temp_t);
        CAROM::Vector temp_y(2,false);
        temp_y(0) = yq1[i];
        temp_y(1) = yq2[i];
        reference_snapshots.push_back(temp_y);
    }

    CAROM::PCHIPInterpolator interp;

    interp.interpolate(times,snapshots,out_times,out_snapshots);

    for(int i = 0; i < out_snapshots.size(); ++i)
    {
        if (abs(reference_snapshots[i](0)-(*out_snapshots[i])(0)) > 10.*FLT_EPSILON)
        {
            SUCCESS = false;
            break;
        }
        else if (abs(reference_snapshots[i](1)-(*out_snapshots[i])(
                         1)) > 10.*FLT_EPSILON)
        {
            SUCCESS = false;
            break;
        }
    }

    out_snapshots.clear();
    out_times.clear();

    interp.interpolate(times,snapshots,26,out_times,out_snapshots);

    for(int i = 0; i < out_snapshots.size(); ++i)
    {
        if( abs(reference_snapshots[i](0)-(*out_snapshots[i])(
                    0)) > 10.*FLT_EPSILON)
        {
            SUCCESS = false;
            break;
        }
        else if(abs(reference_snapshots[i](1)-(*out_snapshots[i])(1)) > 10.*FLT_EPSILON)
        {
            SUCCESS = false;
            break;
        }
    }
    EXPECT_TRUE(SUCCESS);
}


int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
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
