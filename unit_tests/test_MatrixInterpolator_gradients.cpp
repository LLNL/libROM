/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

#include "algo/manifold_interp/MatrixInterpolator.h"
#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include <cfloat>
#include <cmath>
#include <iostream>

#ifndef _WIN32
#include <sys/stat.h>  // mkdir
#else
#include <direct.h>    // _mkdir
#define mkdir(dir, mode) _mkdir(dir)
#endif

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
std::vector<std::shared_ptr<CAROM::Matrix>> createGradient(CAROM::Vector input)
{
    std::vector<double> grad00 = dfunc00(input.getData());
    std::vector<double> grad01 = dfunc01(input.getData());
    std::vector<double> grad10 = dfunc10(input.getData());
    std::vector<double> grad11 = dfunc11(input.getData());

    std::shared_ptr<CAROM::Matrix> gradient0(new Matrix(2,2,false));
    gradient0->item(0,0) = grad00[0]; gradient0->item(0,1) = grad01[0];
    gradient0->item(1,0) = grad10[0]; gradient0->item(1,1) = grad11[0];

    std::shared_ptr<CAROM::Matrix> gradient1(new Matrix(2,2,false));
    gradient1->item(0,0) = grad00[1]; gradient1->item(0,1) = grad01[1];
    gradient1->item(1,0) = grad10[1]; gradient1->item(1,1) = grad11[1];

    std::vector<std::shared_ptr<CAROM::Matrix>> gradient{gradient0, gradient1};

    return gradient;
}

void printMatrix(std::shared_ptr<CAROM::Matrix> input, int rows, int cols, string name)
{
    std::cout << name << "=" << std::endl;
    for(int i = 0; i < rows; ++i)
    {
        for(int j = 0; j < cols; ++j)
        {
            std::cout << input->item(i,j) << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << "-------------------" << std::endl;
}

int main(int argc, char *argv[])
{
    // CAROM::Vector input1(2, false);
    // input1(0) = 0.0; input1(1) = 0.0;
    // CAROM::Vector input2(2, false);
    // input2(0) = 0.0; input2(1) = 0.1;
    // CAROM::Vector input3(2, false);
    // input3(0) = 0.1; input3(1) = 0.0;
    // CAROM::Vector input4(2, false);
    // input4(0) = 0.1; input4(1) = 0.1;
    CAROM::Vector input1(2, false);
    input1(0) = 0.04; input1(1) = 0.04;
    CAROM::Vector input2(2, false);
    input2(0) = 0.04; input2(1) = 0.06;
    CAROM::Vector input3(2, false);
    input3(0) = 0.06; input3(1) = 0.04;
    CAROM::Vector input4(2, false);
    input4(0) = 0.06; input4(1) = 0.06;

    CAROM::Vector target(2,false);
    target(0) = 0.05; target(1) = 0.05;

    //Interpolate a 2x2 matrix.  We should have 4 different functions.
    
    std::shared_ptr<CAROM::Matrix> matrix1 = createMatrix(input1);
    std::shared_ptr<CAROM::Matrix> matrix2 = createMatrix(input2);
    std::shared_ptr<CAROM::Matrix> matrix3 = createMatrix(input3);
    std::shared_ptr<CAROM::Matrix> matrix4 = createMatrix(input4);

    std::shared_ptr<CAROM::Matrix> soln = createMatrix(target);
    

    std::vector<std::shared_ptr<CAROM::Matrix>> true_gradient = createGradient(target);
    
    std::shared_ptr<CAROM::Matrix> eye(new Matrix(2,2,false));
    eye->item(0,0) = 1.0; eye->item(0,1) = 0.0;
    eye->item(0,1) = 0.0; eye->item(1,1) = 1.0;

    std::vector<CAROM::Vector> inputs{input1, input2, input3, input4};
    std::vector<std::shared_ptr<CAROM::Matrix>> rotations{eye,eye,eye,eye};
    std::vector<std::shared_ptr<CAROM::Matrix>> matrices{matrix1, matrix2, matrix3, matrix4};
    
    CAROM::MatrixInterpolator interpolator(inputs, rotations, matrices, 0, "B", "G", "LS", 0.9, true);

    std::shared_ptr<CAROM::Matrix> A = interpolator.interpolate(target);
    

    std::vector<std::shared_ptr<CAROM::Matrix>> gradient = interpolator.getGradient();

    // double epsilon = std::sqrt(std::numeric_limits<float>::denorm_min());
    double epsilon = 1.0e-7;
    CAROM::Vector perturb0(2, false);
    perturb0(0) = epsilon; perturb0(1) = 0.0;
    perturb0 += target;
    CAROM::Vector perturb1(2, false);
    perturb1(0) = 0.0; perturb1(1) = epsilon;
    perturb1 += target;

    std::shared_ptr<CAROM::Matrix> perturbation0 = createMatrix(perturb0);
    std::shared_ptr<CAROM::Matrix> perturbation1 = createMatrix(perturb1);

    std::shared_ptr<CAROM::Matrix> fd_grad0(new Matrix(perturbation0->getData(), 2,2,false));
    *fd_grad0 -= *soln;
    
    std::shared_ptr<CAROM::Matrix> fd_grad1(new Matrix(perturbation1->getData(), 2,2,false));
    *fd_grad1 -= *soln;

    for(int i = 0; i < 2; ++i)
    {
        for(int j = 0; j < 2; ++j)
        {
            fd_grad0->item(i,j) *= 1./epsilon; 
            fd_grad1->item(i,j) *= 1./epsilon;
        }
    }






    

    printMatrix(soln, 2,2, "true solution");
    printMatrix(A,2,2,"interpolated Matrix");

    printMatrix(true_gradient[0],2,2,"true_grad0");
    printMatrix(gradient[0],2,2,"grad0");
    printMatrix(fd_grad0, 2,2, "finite difference grad0");

    printMatrix(true_gradient[1],2,2,"true_grad1");
    printMatrix(gradient[1],2,2,"grad1");
    printMatrix(fd_grad1, 2,2, "finite difference grad1");
    
}