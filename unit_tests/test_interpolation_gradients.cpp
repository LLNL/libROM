/******************************************************************************
 *
 * Copyright (c) 2013-2025, Lawrence Livermore National Security, LLC
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
#include "algo/manifold_interp/MatrixInterpolator.h"
#include "algo/manifold_interp/VectorInterpolator.h"
#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include <cfloat>

using namespace std;
using namespace CAROM;

double func00(double* input)
{
    return 4*input[0] + 3*input[1];
}
std::vector<double> dfunc00(double* input)
{
    std::vector<double> output{4.0, 3.0};
    return output;
}
double func01(double* input)
{
    return 2*input[0]*input[1] + 2*input[0]*input[0];
}
std::vector<double> dfunc01(double* input)
{
    std::vector<double> output{2*input[1]+4*input[0], 2*input[0]};
    return output;
}
double func10(double* input)
{
    return input[1]*input[1] + 2.0;
}
std::vector<double> dfunc10(double* input)
{
    std::vector<double> output{0.0, 2*input[1]};
    return output;
}
double func11(double* input)
{
    return std::cos(3*input[0]);
}
std::vector<double> dfunc11(double* input)
{
    std::vector<double> output{-3*std::sin(3*input[0]), 0.0};
    return output;
}
std::shared_ptr<CAROM::Matrix> createMatrix(CAROM::Vector input)
{
    std::shared_ptr<CAROM::Matrix> output(new Matrix(2,2,false));
    output->item(0,0) = func00(input.getData());
    output->item(0,1) = func01(input.getData());
    output->item(1,0) = func10(input.getData());
    output->item(1,1) = func11(input.getData());
    return output;
}

std::shared_ptr<CAROM::Vector> createVector(CAROM::Vector input)
{
    std::shared_ptr<CAROM::Vector> output(new Vector(2,false));
    output->item(0) = func00(input.getData());
    output->item(1) = func01(input.getData());
    return output;
}

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
    SUCCEED();
}

TEST(InterpolationTest,GaussianMatrix)
{
    bool SUCCESS = true;
    CAROM::Vector input1(2, false);
    input1(0) = 0.04;
    input1(1) = 0.04;
    CAROM::Vector input2(2, false);
    input2(0) = 0.04;
    input2(1) = 0.06;
    CAROM::Vector input3(2, false);
    input3(0) = 0.06;
    input3(1) = 0.04;
    CAROM::Vector input4(2, false);
    input4(0) = 0.06;
    input4(1) = 0.06;

    CAROM::Vector target(2,false);
    target(0) = 0.05;
    target(1) = 0.05;

    //Interpolate a 2x2 matrix.  We should have 4 different functions.

    std::shared_ptr<CAROM::Matrix> matrix1 = createMatrix(input1);
    std::shared_ptr<CAROM::Matrix> matrix2 = createMatrix(input2);
    std::shared_ptr<CAROM::Matrix> matrix3 = createMatrix(input3);
    std::shared_ptr<CAROM::Matrix> matrix4 = createMatrix(input4);

    std::shared_ptr<CAROM::Matrix> soln = createMatrix(target);

    std::shared_ptr<CAROM::Matrix> eye(new Matrix(2,2,false));
    eye->item(0,0) = 1.0;
    eye->item(0,1) = 0.0;
    eye->item(0,1) = 0.0;
    eye->item(1,1) = 1.0;

    std::vector<CAROM::Vector> inputs{input1, input2, input3, input4};
    std::vector<std::shared_ptr<CAROM::Matrix>> rotations{eye,eye,eye,eye};
    std::vector<std::shared_ptr<CAROM::Matrix>> matrices{matrix1, matrix2, matrix3, matrix4};

    CAROM::MatrixInterpolator interpolator(inputs, rotations, matrices, 0, "B", "G",
                                           "LS", 0.9, true);

    std::shared_ptr<CAROM::Matrix> A = interpolator.interpolate(target);
    std::vector<std::shared_ptr<CAROM::Matrix>> gradient =
                interpolator.getGradient();

    double epsilon = 1.0e-7;
    CAROM::Vector perturb0(2, false);
    perturb0(0) = epsilon;
    perturb0(1) = 0.0;
    perturb0 += target;
    CAROM::Vector perturb1(2, false);
    perturb1(0) = 0.0;
    perturb1(1) = epsilon;
    perturb1 += target;

    std::shared_ptr<CAROM::Matrix> perturbation0 = createMatrix(perturb0);
    std::shared_ptr<CAROM::Matrix> perturbation1 = createMatrix(perturb1);

    std::shared_ptr<CAROM::Matrix> interp_perturbation0 = interpolator.interpolate(
                perturb0);
    std::shared_ptr<CAROM::Matrix> interp_perturbation1 = interpolator.interpolate(
                perturb1);

    std::shared_ptr<CAROM::Matrix> fd_grad_interp0(new Matrix(
                *interp_perturbation0));
    *fd_grad_interp0 -= *A;

    std::shared_ptr<CAROM::Matrix> fd_grad_interp1(new Matrix(
                *interp_perturbation1));
    *fd_grad_interp1 -= *A;

    for(int i = 0; i < 2; ++i)
    {
        for(int j = 0; j < 2; ++j)
        {
            fd_grad_interp0->item(i,j) /= epsilon;
            fd_grad_interp1->item(i,j) /= epsilon;
        }
    }

    for(int i = 0; i < 2; ++i)
    {
        for(int j = 0; j < 2; ++j)
        {
            double fd = fd_grad_interp0->item(i,j);
            double my = gradient[0]->item(i,j);
            if(abs(fd - my) > 1.e-5)
            {
                SUCCESS = false;
                std::cout << "fd = " << fd << ", my = " << my << "abs(fd-my) = " << abs(
                              fd-my) << std::endl;
                break;
            }
            fd = fd_grad_interp1->item(i,j);
            my = gradient[1]->item(i,j);
            if(abs(fd - my) > 1.e-5)
            {
                SUCCESS = false;
                std::cout << "fd = " << fd << ", my = " << my << "abs(fd-my) = " << abs(
                              fd-my) << std::endl;
                break;
            }
        }
        if(!SUCCESS) break;
    }
    EXPECT_TRUE(SUCCESS);
}
TEST(InterpolationTest,InverseQuadraticMatrix)
{
    bool SUCCESS = true;
    CAROM::Vector input1(2, false);
    input1(0) = 0.04;
    input1(1) = 0.04;
    CAROM::Vector input2(2, false);
    input2(0) = 0.04;
    input2(1) = 0.06;
    CAROM::Vector input3(2, false);
    input3(0) = 0.06;
    input3(1) = 0.04;
    CAROM::Vector input4(2, false);
    input4(0) = 0.06;
    input4(1) = 0.06;

    CAROM::Vector target(2,false);
    target(0) = 0.05;
    target(1) = 0.05;

    //Interpolate a 2x2 matrix.  We should have 4 different functions.

    std::shared_ptr<CAROM::Matrix> matrix1 = createMatrix(input1);
    std::shared_ptr<CAROM::Matrix> matrix2 = createMatrix(input2);
    std::shared_ptr<CAROM::Matrix> matrix3 = createMatrix(input3);
    std::shared_ptr<CAROM::Matrix> matrix4 = createMatrix(input4);

    std::shared_ptr<CAROM::Matrix> soln = createMatrix(target);

    std::shared_ptr<CAROM::Matrix> eye(new Matrix(2,2,false));
    eye->item(0,0) = 1.0;
    eye->item(0,1) = 0.0;
    eye->item(0,1) = 0.0;
    eye->item(1,1) = 1.0;

    std::vector<CAROM::Vector> inputs{input1, input2, input3, input4};
    std::vector<std::shared_ptr<CAROM::Matrix>> rotations{eye,eye,eye,eye};
    std::vector<std::shared_ptr<CAROM::Matrix>> matrices{matrix1, matrix2, matrix3, matrix4};

    CAROM::MatrixInterpolator interpolator(inputs, rotations, matrices, 0, "B",
                                           "IQ", "LS", 0.9, true);

    std::shared_ptr<CAROM::Matrix> A = interpolator.interpolate(target);
    std::vector<std::shared_ptr<CAROM::Matrix>> gradient =
                interpolator.getGradient();

    double epsilon = 1.0e-7;
    CAROM::Vector perturb0(2, false);
    perturb0(0) = epsilon;
    perturb0(1) = 0.0;
    perturb0 += target;
    CAROM::Vector perturb1(2, false);
    perturb1(0) = 0.0;
    perturb1(1) = epsilon;
    perturb1 += target;

    std::shared_ptr<CAROM::Matrix> perturbation0 = createMatrix(perturb0);
    std::shared_ptr<CAROM::Matrix> perturbation1 = createMatrix(perturb1);

    std::shared_ptr<CAROM::Matrix> interp_perturbation0 = interpolator.interpolate(
                perturb0);
    std::shared_ptr<CAROM::Matrix> interp_perturbation1 = interpolator.interpolate(
                perturb1);

    std::shared_ptr<CAROM::Matrix> fd_grad_interp0(new Matrix(
                *interp_perturbation0));
    *fd_grad_interp0 -= *A;

    std::shared_ptr<CAROM::Matrix> fd_grad_interp1(new Matrix(
                *interp_perturbation1));
    *fd_grad_interp1 -= *A;

    for(int i = 0; i < 2; ++i)
    {
        for(int j = 0; j < 2; ++j)
        {
            fd_grad_interp0->item(i,j) /= epsilon;
            fd_grad_interp1->item(i,j) /= epsilon;
        }
    }

    for(int i = 0; i < 2; ++i)
    {
        for(int j = 0; j < 2; ++j)
        {
            double fd = fd_grad_interp0->item(i,j);
            double my = gradient[0]->item(i,j);
            if(abs(fd - my) > 1.e-5)
            {
                SUCCESS = false;
                break;
            }
            fd = fd_grad_interp0->item(i,j);
            my = gradient[0]->item(i,j);
            if(abs(fd - my) > 1.e-5)
            {
                SUCCESS = false;
                break;
            }
        }
        if(!SUCCESS) break;
    }
    EXPECT_TRUE(SUCCESS);
}

TEST(InterpolationTest,InverseMultiQuadraticMatrix)
{
    bool SUCCESS = true;
    CAROM::Vector input1(2, false);
    input1(0) = 0.04;
    input1(1) = 0.04;
    CAROM::Vector input2(2, false);
    input2(0) = 0.04;
    input2(1) = 0.06;
    CAROM::Vector input3(2, false);
    input3(0) = 0.06;
    input3(1) = 0.04;
    CAROM::Vector input4(2, false);
    input4(0) = 0.06;
    input4(1) = 0.06;

    CAROM::Vector target(2,false);
    target(0) = 0.05;
    target(1) = 0.05;

    //Interpolate a 2x2 matrix.  We should have 4 different functions.

    std::shared_ptr<CAROM::Matrix> matrix1 = createMatrix(input1);
    std::shared_ptr<CAROM::Matrix> matrix2 = createMatrix(input2);
    std::shared_ptr<CAROM::Matrix> matrix3 = createMatrix(input3);
    std::shared_ptr<CAROM::Matrix> matrix4 = createMatrix(input4);

    std::shared_ptr<CAROM::Matrix> soln = createMatrix(target);

    std::shared_ptr<CAROM::Matrix> eye(new Matrix(2,2,false));
    eye->item(0,0) = 1.0;
    eye->item(0,1) = 0.0;
    eye->item(0,1) = 0.0;
    eye->item(1,1) = 1.0;

    std::vector<CAROM::Vector> inputs{input1, input2, input3, input4};
    std::vector<std::shared_ptr<CAROM::Matrix>> rotations{eye,eye,eye,eye};
    std::vector<std::shared_ptr<CAROM::Matrix>> matrices{matrix1, matrix2, matrix3, matrix4};

    CAROM::MatrixInterpolator interpolator(inputs, rotations, matrices, 0, "B",
                                           "IMQ", "LS", 0.9, true);

    std::shared_ptr<CAROM::Matrix> A = interpolator.interpolate(target);
    std::vector<std::shared_ptr<CAROM::Matrix>> gradient =
                interpolator.getGradient();

    double epsilon = 1.0e-7;
    CAROM::Vector perturb0(2, false);
    perturb0(0) = epsilon;
    perturb0(1) = 0.0;
    perturb0 += target;
    CAROM::Vector perturb1(2, false);
    perturb1(0) = 0.0;
    perturb1(1) = epsilon;
    perturb1 += target;

    std::shared_ptr<CAROM::Matrix> perturbation0 = createMatrix(perturb0);
    std::shared_ptr<CAROM::Matrix> perturbation1 = createMatrix(perturb1);

    std::shared_ptr<CAROM::Matrix> interp_perturbation0 = interpolator.interpolate(
                perturb0);
    std::shared_ptr<CAROM::Matrix> interp_perturbation1 = interpolator.interpolate(
                perturb1);

    std::shared_ptr<CAROM::Matrix> fd_grad_interp0(new Matrix(
                *interp_perturbation0));
    *fd_grad_interp0 -= *A;

    std::shared_ptr<CAROM::Matrix> fd_grad_interp1(new Matrix(
                *interp_perturbation1));
    *fd_grad_interp1 -= *A;

    for(int i = 0; i < 2; ++i)
    {
        for(int j = 0; j < 2; ++j)
        {
            fd_grad_interp0->item(i,j) /= epsilon;
            fd_grad_interp1->item(i,j) /= epsilon;
        }
    }

    for(int i = 0; i < 2; ++i)
    {
        for(int j = 0; j < 2; ++j)
        {
            double fd = fd_grad_interp0->item(i,j);
            double my = gradient[0]->item(i,j);
            if(abs(fd - my) > 1.e-5)
            {
                SUCCESS = false;
                break;
            }
            fd = fd_grad_interp0->item(i,j);
            my = gradient[0]->item(i,j);
            if(abs(fd - my) > 1.e-5)
            {
                SUCCESS = false;
                break;
            }
        }
        if(!SUCCESS) break;
    }
    EXPECT_TRUE(SUCCESS);
}

TEST(InterpolationTest,Vector)
{
    bool SUCCESS = true;
    CAROM::Vector input1(2, false);
    input1(0) = 0.04;
    input1(1) = 0.04;
    CAROM::Vector input2(2, false);
    input2(0) = 0.04;
    input2(1) = 0.06;
    CAROM::Vector input3(2, false);
    input3(0) = 0.06;
    input3(1) = 0.04;
    CAROM::Vector input4(2, false);
    input4(0) = 0.06;
    input4(1) = 0.06;

    CAROM::Vector target(2,false);
    target(0) = 0.05;
    target(1) = 0.05;


    std::shared_ptr<CAROM::Vector> matrix1 = createVector(input1);
    std::shared_ptr<CAROM::Vector> matrix2 = createVector(input2);
    std::shared_ptr<CAROM::Vector> matrix3 = createVector(input3);
    std::shared_ptr<CAROM::Vector> matrix4 = createVector(input4);

    std::shared_ptr<CAROM::Vector> soln = createVector(target);

    std::shared_ptr<CAROM::Matrix> eye(new Matrix(2,2,false));
    eye->item(0,0) = 1.0;
    eye->item(0,1) = 0.0;
    eye->item(0,1) = 0.0;
    eye->item(1,1) = 1.0;

    std::vector<CAROM::Vector> inputs{input1, input2, input3, input4};
    std::vector<std::shared_ptr<CAROM::Matrix>> rotations{eye,eye,eye,eye};
    std::vector<std::shared_ptr<CAROM::Vector>> vectors{matrix1, matrix2, matrix3, matrix4};

    CAROM::VectorInterpolator interpolator(inputs, rotations, vectors, 0, "G", "LS",
                                           0.9, true);

    std::shared_ptr<CAROM::Vector> v = interpolator.interpolate(target);
    std::vector<std::shared_ptr<CAROM::Vector>> gradient =
                interpolator.getGradient();

    double epsilon = 1.0e-7;
    CAROM::Vector perturb0(2, false);
    perturb0(0) = epsilon;
    perturb0(1) = 0.0;
    perturb0 += target;
    CAROM::Vector perturb1(2, false);
    perturb1(0) = 0.0;
    perturb1(1) = epsilon;
    perturb1 += target;

    std::shared_ptr<CAROM::Vector> perturbation0 = createVector(perturb0);
    std::shared_ptr<CAROM::Vector> perturbation1 = createVector(perturb1);

    std::shared_ptr<CAROM::Vector> interp_perturbation0 = interpolator.interpolate(
                perturb0);
    std::shared_ptr<CAROM::Vector> interp_perturbation1 = interpolator.interpolate(
                perturb1);

    std::shared_ptr<CAROM::Vector> fd_grad_interp0(new Vector(
                *interp_perturbation0));
    *fd_grad_interp0 -= *v;

    std::shared_ptr<CAROM::Vector> fd_grad_interp1(new Vector(
                *interp_perturbation1));
    *fd_grad_interp1 -= *v;

    for(int i = 0; i < 2; ++i)
    {
        fd_grad_interp0->item(i) /= epsilon;
        fd_grad_interp1->item(i) /= epsilon;
    }

    for(int i = 0; i < 2; ++i)
    {
        double fd = fd_grad_interp0->item(i);
        double my = gradient[0]->item(i);
        if(abs(fd - my) > 1.e-5)
        {
            SUCCESS = false;
            break;
        }
        fd = fd_grad_interp0->item(i);
        my = gradient[0]->item(i);
        if(abs(fd - my) > 1.e-5)
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
